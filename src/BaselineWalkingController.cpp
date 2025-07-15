#include <sys/syscall.h>

#include <mc_tasks/CoMTask.h>
#include <mc_tasks/FirstOrderImpedanceTask.h>
#include <mc_tasks/MetaTaskLoader.h>
#include <mc_tasks/OrientationTask.h>

#include <BaselineWalkingController/BaselineWalkingController.h>
#include <BaselineWalkingController/CentroidalManager.h>
#include <BaselineWalkingController/FootManager.h>
#include <BaselineWalkingController/centroidal/CentroidalManagerDdpZmp.h>
#include <BaselineWalkingController/centroidal/CentroidalManagerFootGuidedControl.h>
#include <BaselineWalkingController/centroidal/CentroidalManagerIntrinsicallyStableMpc.h>
#include <BaselineWalkingController/centroidal/CentroidalManagerPreviewControlZmp.h>

using namespace BWC;

BaselineWalkingController::BaselineWalkingController(mc_rbdyn::RobotModulePtr rm,
                                                     double dt,
                                                     const mc_rtc::Configuration & _config,
                                                     bool allowEmptyManager)
: mc_control::fsm::Controller(rm, dt, _config)
{
  // ロボット固有の設定を取得
  auto rconfig = config()("robots")(robot().module().name);
  if(rconfig.empty())
  {
    mc_rtc::log::error_and_throw("[BaselineWalkingController] {} section is empty, please provide a configuration",
                                 robot().module().name);
  }
  // 設定を読み込み
  config().load(rconfig);

  // 追加の上書き設定を読み込み
  auto overwriteConfigList = config()("OverwriteConfigList", mc_rtc::Configuration());
  auto overwriteConfigKeys = config()("OverwriteConfigKeys", std::vector<std::string>{});
  for(const auto & overwriteConfigKey : overwriteConfigKeys)
  {
    if(!overwriteConfigList.has(overwriteConfigKey))
    {
      mc_rtc::log::error_and_throw(
          "[BaselineWalkingController] {} in OverwriteConfigKeys but not in OverwriteConfigList", overwriteConfigKey);
    }
    config().load(overwriteConfigList(overwriteConfigKey));
  }

  config()("controllerName", name_);

  // CoMタスクのセットアップ
  if(config().has("CoMTask"))
  {
    comTask_ = mc_tasks::MetaTaskLoader::load<mc_tasks::CoMTask>(solver(), config()("CoMTask"));
    comTask_->name("CoMTask");
  }
  else
  {
    mc_rtc::log::warning("[BaselineWalkingController] CoMTask configuration is missing.");
  }

  // 基底姿勢タスクのセットアップ
  if(config().has("BaseOrientationTask"))
  {
    baseOriTask_ = mc_tasks::MetaTaskLoader::load<mc_tasks::OrientationTask>(solver(), config()("BaseOrientationTask"));
    baseOriTask_->name("BaseOriTask");
  }
  else
  {
    mc_rtc::log::warning("[BaselineWalkingController] BaseOrientationTask configuration is missing.");
  }

  // 足タスクのセットアップ
  if(config().has("FootTaskList"))
  {
    for(const auto & footTaskConfig : config()("FootTaskList"))
    {
      Foot foot = strToFoot(footTaskConfig("foot"));
      footTasks_.emplace(
          foot, mc_tasks::MetaTaskLoader::load<mc_tasks::force::FirstOrderImpedanceTask>(solver(), footTaskConfig));
      footTasks_.at(foot)->name("FootTask_" + std::to_string(foot));
    }
  }
  else
  {
    mc_rtc::log::warning("[BaselineWalkingController] FootTaskList configuration is missing.");
  }

  // FootManagerのセットアップ
  if(config().has("FootManager"))
  {
    footManager_ = std::make_shared<FootManager>(this, config()("FootManager"));
  }
  else
  {
    mc_rtc::log::warning("[BaselineWalkingController] FootManager configuration is missing.");
  }

  // CentroidalManagerのセットアップ
  if(config().has("CentroidalManager"))
  {
    std::string centroidalManagerMethod = config()("CentroidalManager")("method", std::string(""));
    if(centroidalManagerMethod == "PreviewControlZmp")
    {
      centroidalManager_ = std::make_shared<CentroidalManagerPreviewControlZmp>(this, config()("CentroidalManager"));
    }
    else if(centroidalManagerMethod == "DdpZmp")
    {
      centroidalManager_ = std::make_shared<CentroidalManagerDdpZmp>(this, config()("CentroidalManager"));
    }
    else if(centroidalManagerMethod == "FootGuidedControl")
    {
      centroidalManager_ = std::make_shared<CentroidalManagerFootGuidedControl>(this, config()("CentroidalManager"));
    }
    else if(centroidalManagerMethod == "IntrinsicallyStableMpc")
    {
      centroidalManager_ =
          std::make_shared<CentroidalManagerIntrinsicallyStableMpc>(this, config()("CentroidalManager"));
    }
    else
    {
      if(!allowEmptyManager)
      {
        mc_rtc::log::error_and_throw("[BaselineWalkingController] Invalid centroidalManagerMethod: {}.",
                                     centroidalManagerMethod);
      }
    }
  }
  else
  {
    mc_rtc::log::warning("[BaselineWalkingController] CentroidalManager configuration is missing.");
  }

  // デフォルトアンカーのセット
  setDefaultAnchor();

  mc_rtc::log::success("[BaselineWalkingController] Constructed.");
}

void BaselineWalkingController::reset(const mc_control::ControllerResetData & resetData)
{
  mc_control::fsm::Controller::reset(resetData);

  enableManagerUpdate_ = false;

  // スレッドIDを取得し、優先度変更方法を表示
  long tid = static_cast<long>(syscall(SYS_gettid));
  mc_rtc::log::info("[BaselineWalkingController] TID is {}. Run the following command to set high priority:\n  sudo "
                    "renice -n -20 -p {}",
                    tid, tid);
  mc_rtc::log::info("[BaselineWalkingController] You can check the current priority by the following command:\n  ps -p "
                    "`pgrep choreonoid` -o pid,tid,args,ni,pri,wchan m");

  mc_rtc::log::success("[BaselineWalkingController] Reset.");
}

bool BaselineWalkingController::run()
{
  t_ += dt();

  if(enableManagerUpdate_)
  {
    // 各マネージャを更新
    footManager_->update();
    centroidalManager_->update();
  }

  return mc_control::fsm::Controller::run();
}

void BaselineWalkingController::stop()
{
  // 各タスクを削除
  solver().removeTask(comTask_);
  solver().removeTask(baseOriTask_);
  for(const auto & foot : Feet::Both)
  {
    solver().removeTask(footTasks_.at(foot));
  }

  // 各マネージャを停止
  footManager_->stop();
  centroidalManager_->stop();

  // デフォルトアンカーに戻す
  setDefaultAnchor();

  mc_control::fsm::Controller::stop();
}

void BaselineWalkingController::setDefaultAnchor()
{
  // 左右足の中点をアンカーとする
  std::string anchorName = "KinematicAnchorFrame::" + robot().name();
  if(datastore().has(anchorName))
  {
    datastore().remove(anchorName);
  }
  datastore().make_call(anchorName, [this](const mc_rbdyn::Robot & robot) {
    return sva::interpolate(robot.surfacePose(footManager_->surfaceName(Foot::Left)),
                            robot.surfacePose(footManager_->surfaceName(Foot::Right)), 0.5);
  });
}
