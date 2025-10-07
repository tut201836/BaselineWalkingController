#include <RBDyn/Momentum.h> //add
#include <fstream>

#include <mc_rtc/gui/Button.h>
#include <mc_rtc/gui/Checkbox.h>
#include <mc_rtc/gui/ComboInput.h>
#include <mc_rtc/gui/Label.h>
#include <mc_rtc/gui/NumberInput.h>
#include <mc_rtc/gui/plot.h>
#include <mc_tasks/CoMTask.h>
#include <mc_tasks/FirstOrderImpedanceTask.h>

#include <CCC/Constants.h>
#include <ForceColl/WrenchDistribution.h>

#include <BaselineWalkingController/BaselineWalkingController.h>
#include <BaselineWalkingController/CentroidalManager.h>
#include <BaselineWalkingController/FootManager.h>

using namespace BWC;

bool isDangerFloor = false;

void CentroidalManager::DcmEstimatorConfiguration::load(const mc_rtc::Configuration & mcRtcConfig)
{
  mcRtcConfig("enableDcmEstimator", enableDcmEstimator);
  mcRtcConfig("dcmCorrectionMode", dcmCorrectionMode);
  const std::vector<std::string> dcmCorrectionModeList = {"Bias", "Filter", "None"};
  if(std::find(dcmCorrectionModeList.begin(), dcmCorrectionModeList.end(), dcmCorrectionMode)
     == dcmCorrectionModeList.end())
  {
    mc_rtc::log::error_and_throw("[CentroidalManager] Invalid dcmCorrectionMode: {}", dcmCorrectionMode);
  }
  mcRtcConfig("biasDriftPerSecondStd", biasDriftPerSecondStd);
  mcRtcConfig("dcmMeasureErrorStd", dcmMeasureErrorStd);
  mcRtcConfig("zmpMeasureErrorStd", zmpMeasureErrorStd);
  mcRtcConfig("biasLimit", biasLimit);
  mcRtcConfig("initDcmUncertainty", initDcmUncertainty);
  mcRtcConfig("initBiasUncertainty", initBiasUncertainty);
}

void CentroidalManager::Configuration::load(const mc_rtc::Configuration & mcRtcConfig)
{
  mcRtcConfig("name", name);
  mcRtcConfig("method", method);
  mcRtcConfig("useActualStateForMpc", useActualStateForMpc);
  mcRtcConfig("enableZmpFeedback", enableZmpFeedback);
  mcRtcConfig("enableComZFeedback", enableComZFeedback);
  mcRtcConfig("dcmGainP", dcmGainP);
  mcRtcConfig("zmpVelGain", zmpVelGain);
  mcRtcConfig("comZGainP", comZGainP);
  mcRtcConfig("comZGainD", comZGainD);
  mcRtcConfig("refComZ", refComZ);
  mcRtcConfig("useTargetPoseForControlRobotAnchorFrame", useTargetPoseForControlRobotAnchorFrame);
  mcRtcConfig("useActualComForWrenchDist", useActualComForWrenchDist);
  mcRtcConfig("actualComOffset", actualComOffset);
  mcRtcConfig("wrenchDistConfig", wrenchDistConfig);
  mcRtcConfig("floorSpringK", floorSpringK);
  mcRtcConfig("floorDampingD", floorDampingD);
  mcRtcConfig("sinkOnThreshold", sinkOnThreshold);   // default 0.010
  mcRtcConfig("sinkOffThreshold", sinkOffThreshold); // default 0.006
  mcRtcConfig("sinkAlpha", sinkAlpha);               // LPF alpha (0..1), default 0.25
  mcRtcConfig("sinkAccelScale", sinkAccelScale);     // 既存の 0.001 を置換
  dcmEstimatorConfig.load(mcRtcConfig("dcmEstimator", mc_rtc::Configuration()));
}

CentroidalManager::CentroidalManager(BaselineWalkingController * ctlPtr, const mc_rtc::Configuration & // mcRtcConfig
                                     )
: ctlPtr_(ctlPtr), refComZFunc_(std::make_shared<TrajColl::CubicInterpolator<double>>())
{
}

void CentroidalManager::reset()
{
  robotMass_ = ctl().robot().mass();

  refComZFunc_->clearPoints();
  refComZFunc_->appendPoint(std::make_pair(ctl().t(), config().refComZ));
  refComZFunc_->appendPoint(std::make_pair(interpMaxTime_, config().refComZ));
  refComZFunc_->calcCoeff();

  const DcmEstimatorConfiguration & dcmEstimatorConfig = config().dcmEstimatorConfig;
  double nominalOmega = std::sqrt(mc_rtc::constants::GRAVITY / config().refComZ);
  dcmEstimator_ = std::make_shared<stateObservation::LipmDcmEstimator>(
      ctl().dt(), nominalOmega, dcmEstimatorConfig.biasDriftPerSecondStd, dcmEstimatorConfig.dcmMeasureErrorStd,
      dcmEstimatorConfig.zmpMeasureErrorStd, dcmEstimatorConfig.biasLimit, Eigen::Vector2d::Zero(),
      Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero(), dcmEstimatorConfig.initDcmUncertainty,
      dcmEstimatorConfig.initBiasUncertainty);
  requireDcmEstimatorReset_ = true;
}

void CentroidalManager::update()
{
  // Set MPC state
  if(config().useActualStateForMpc)
  {
    mpcCom_ = actualCom();
    mpcComVel_ = ctl().realRobot().comVelocity();
  }
  else
  {
    // Task targets are the planned state in the previous step
    mpcCom_ = ctl().comTask_->com();
    mpcComVel_ = ctl().comTask_->refVel();
  }
  refZmp_ = ctl().footManager_->calcRefZmp(ctl().t());

  // Run MPC
  runMpc();

  // Calculate measured ZMP
  {
    std::unordered_map<Foot, sva::ForceVecd> sensorWrenchList;
    for(const auto & foot : ctl().footManager_->getCurrentContactFeet())
    {
      const auto & surfaceName = ctl().footManager_->surfaceName(foot);
      const auto & sensorName = ctl().robot().indirectSurfaceForceSensor(surfaceName).name();
      const auto & sensor = ctl().robot().forceSensor(sensorName);
      const auto & sensorWrench = sensor.worldWrenchWithoutGravity(ctl().robot());
      sensorWrenchList.emplace(foot, sensorWrench);
    }
    measuredZMP_ = calcZmp(sensorWrenchList, refZmp_.z());
  }

  // Calculate target wrench
  {
    controlZmp_ = plannedZmp_;
    controlForceZ_ = plannedForceZ_;

    // Compensate ZMP delay
    // See equation (10) of https://ieeexplore.ieee.org/abstract/document/6094838
    Eigen::Vector3d refZmpVel = ctl().footManager_->calcRefZmp(ctl().t(), 1);
    controlZmp_.head<2>() += config().zmpVelGain * refZmpVel.head<2>();

    // Apply DCM feedback
    if(config().enableZmpFeedback)
    {
      double omega = std::sqrt(plannedForceZ_ / (robotMass_ * (mpcCom_.z() - refZmp_.z())));
      Eigen::Vector3d plannedDcm = ctl().comTask_->com() + ctl().comTask_->refVel() / omega;
      Eigen::Vector3d actualDcm = actualCom() + ctl().realRobot().comVelocity() / omega;
      Eigen::Vector3d dcmError = actualDcm - plannedDcm;


      // Estimate and compensate for DCM bias
      if(config().dcmEstimatorConfig.enableDcmEstimator)
      {
        if(omega != dcmEstimator_->getLipmNaturalFrequency())
        {
          dcmEstimator_->setLipmNaturalFrequency(omega);
        }

        const Eigen::Matrix3d & baseOrientation = ctl().robot().posW().rotation().transpose();
        if(requireDcmEstimatorReset_)
        {
          dcmEstimator_->resetWithMeasurements(actualDcm.head<2>(), measuredZMP_.head<2>(), baseOrientation, true);
          requireDcmEstimatorReset_ = false;
        }
        else
        {
          dcmEstimator_->setInputs(actualDcm.head<2>(), measuredZMP_.head<2>(), baseOrientation);
        }

        dcmEstimator_->update();

        if(config().dcmEstimatorConfig.dcmCorrectionMode == "Bias")
        {
          dcmError.head<2>() -= dcmEstimator_->getBias();
        }
        else if(config().dcmEstimatorConfig.dcmCorrectionMode == "Filter")
        {
          dcmError.head<2>() = dcmEstimator_->getUnbiasedDCM();
        }
        // else if(config().dcmEstimatorConfig.dcmCorrectionMode == "None")
        // {
        // }
      }
      else
      {
        requireDcmEstimatorReset_ = true;
        dcmEstimator_->setBias(Eigen::Vector2d::Zero());
      }
      controlZmp_.head<2>() += config().dcmGainP * dcmError.head<2>();
    }

    // --- added: roll-rate damper to mitigate lateral sway ---
    {
      // ベースのロール角速度 ωx を取得（IMU/BodySensor 名は環境に合わせて）
      // mc_rtc では realRobot().bodySensor() / robot().bodySensor() のどちらかを使用
      const auto & bs = ctl().realRobot().bodySensor();
      const double omega_roll = bs.angularVelocity().x(); // [rad/s], +は左肩上がり方向
      // 係数（初期目安）
      static const double k_roll_damp = 0.02; // [m / (rad/s)] 0.01〜0.03で調整
      // ロールが左へ倒れる（ωx>0）時はZMPを左（+y）に寄せる ⇒ 符号は -?
      // 右手系定義に依るが、多くのモデルで "ZMPをロール速度と逆向きに" が安定
      controlZmp_.y() += -k_roll_damp * omega_roll;
    }

    // Apply ForceZ feedback
    if(config().enableComZFeedback)
    {
      double plannedComZ = ctl().comTask_->com().z();
      double actualComZ = actualCom().z();
      double plannedComVelZ = ctl().comTask_->refVel().z();
      double actualComVelZ = ctl().realRobot().comVelocity().z();
      controlForceZ_ -=
          config().comZGainP * (actualComZ - plannedComZ) + config().comZGainD * (actualComVelZ - plannedComVelZ);
    }

    // Convert ZMP to wrench and distribute
    contactList_ = ctl().footManager_->calcCurrentContactList();
    wrenchDist_ = std::make_shared<ForceColl::WrenchDistribution>(ForceColl::getContactVecFromMap(contactList_),
                                                                  config().wrenchDistConfig);
    Eigen::Vector3d comForWrenchDist =
        (config().useActualComForWrenchDist ? actualComUnbiased() : ctl().comTask_->com());
    sva::ForceVecd controlWrench;
    controlWrench.force() << controlForceZ_ / (comForWrenchDist.z() - refZmp_.z())
                                 * (comForWrenchDist.head<2>() - controlZmp_.head<2>()),
        controlForceZ_;
    controlWrench.moment().setZero(); // Moment is represented around CoM
    wrenchDist_->run(controlWrench, comForWrenchDist);
  }

  // Set target of tasks
  {
    // Set target of CoM task
    Eigen::Vector3d plannedComAccel = calcPlannedComAccel();
    Eigen::Vector3d nextPlannedCom =
        mpcCom_ + ctl().dt() * mpcComVel_ + 0.5 * std::pow(ctl().dt(), 2) * plannedComAccel;
    Eigen::Vector3d nextPlannedComVel = mpcComVel_ + ctl().dt() * plannedComAccel;

    double footSurfaceDiff = estimateSupportSink();

    // LPF＋ヒステリシス
    config().footSurfaceDiffFilt_ = config().sinkAlpha * footSurfaceDiff
                                   + (1.0 - config().sinkAlpha) * config().footSurfaceDiffFilt_;
    if(!config().sinkActive_ && config().footSurfaceDiffFilt_ > config().sinkOnThreshold)  config().sinkActive_ = true;
    else if(config().sinkActive_ && config().footSurfaceDiffFilt_ < config().sinkOffThreshold) config().sinkActive_ = false;
    isDangerFloor = config().sinkActive_;
        
    // 沈む床を踏んだとき (今は使わない)
    if (config().sinkActive_)
    {
      // mc_rtc::log::warning("footSurfaceDiff: {}", footSurfaceDiff);
      isDangerFloor = true;
      if(isConstantComZ())
      {
        plannedComAccel.z() = calcRefComZ(ctl().t(), 2) + ctl().footManager_->calcRefGroundPosZ(ctl().t(), 2);
        plannedComAccel.z() += config().sinkAccelScale
                         * compliantFloorCorrection(config().footSurfaceDiffFilt_);
      }
      nextPlannedCom = mpcCom_ + ctl().dt() * mpcComVel_ + 0.5 * std::pow(ctl().dt(), 2) * plannedComAccel;
      nextPlannedComVel = mpcComVel_ + ctl().dt() * plannedComAccel;

      const double delta_h = config().footSurfaceDiffFilt_; // [m]（正で沈み）
      const double gamma = 0.7;                             // 反映率（0.5〜0.8推奨）
      if(true) {
        // Foot swing = ctl().footManager_->swingFoot();
        // ctl().footManager_->setNextLandingZOffset(swing, -gamma * delta_h);
      }
    
    }
    // 通常の床
    else
    {
      if(isConstantComZ())
      {
        nextPlannedCom.z() = calcRefComZ(ctl().t()) + ctl().footManager_->calcRefGroundPosZ(ctl().t());
        nextPlannedComVel.z() = calcRefComZ(ctl().t(), 1) + ctl().footManager_->calcRefGroundPosZ(ctl().t(), 1);
        plannedComAccel.z() = calcRefComZ(ctl().t(), 2) + ctl().footManager_->calcRefGroundPosZ(ctl().t(), 2);
      }
    }

    ctl().comTask_->com(nextPlannedCom);
    ctl().comTask_->refVel(nextPlannedComVel);
    ctl().comTask_->refAccel(plannedComAccel);

    // Set target wrench of foot tasks
    const auto & targetWrenchList = ForceColl::calcWrenchList(contactList_, wrenchDist_->resultWrenchRatio_);
    for(const auto & foot : Feet::Both)
    {
      sva::ForceVecd targetWrench = sva::ForceVecd::Zero();
      if(targetWrenchList.count(foot))
      {
        targetWrench = targetWrenchList.at(foot);
      }
      ctl().footTasks_.at(foot)->targetWrenchW(targetWrench);
    }
  }

  // Calculate support region for log
  {
    supportRegion_[0].setConstant(std::numeric_limits<double>::max());
    supportRegion_[1].setConstant(std::numeric_limits<double>::lowest());
    for(const auto & contactKV : contactList_)
    {
      for(const auto & vertexWithRidge : contactKV.second->vertexWithRidgeList_)
      {
        supportRegion_[0] = supportRegion_[0].cwiseMin(vertexWithRidge.vertex.head<2>());
        supportRegion_[1] = supportRegion_[1].cwiseMax(vertexWithRidge.vertex.head<2>());
      }
    }
  }

  // Update force visualization
  {
    ctl().gui()->removeCategory({ctl().name(), config().name, "ForceMarker"});
    wrenchDist_->addToGUI(*ctl().gui(), {ctl().name(), config().name, "ForceMarker"});
  }
}

void CentroidalManager::stop()
{
  removeFromGUI(*ctl().gui());
  removeFromLogger(ctl().logger());
}

void CentroidalManager::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  gui.addElement(
      {ctl().name(), config().name, "Config"}, mc_rtc::gui::Label("method", [this]() { return config().method; }),
      mc_rtc::gui::Checkbox(
          "useActualStateForMpc", [this]() { return config().useActualStateForMpc; },
          [this]() { config().useActualStateForMpc = !config().useActualStateForMpc; }),
      mc_rtc::gui::Checkbox(
          "enableZmpFeedback", [this]() { return config().enableZmpFeedback; },
          [this]() { config().enableZmpFeedback = !config().enableZmpFeedback; }),
      mc_rtc::gui::Checkbox(
          "enableComZFeedback", [this]() { return config().enableComZFeedback; },
          [this]() { config().enableComZFeedback = !config().enableComZFeedback; }),
      mc_rtc::gui::NumberInput(
          "dcmGainP", [this]() { return config().dcmGainP; }, [this](double v) { config().dcmGainP = v; }),
      mc_rtc::gui::NumberInput(
          "zmpVelGain", [this]() { return config().zmpVelGain; }, [this](double v) { config().zmpVelGain = v; }),
      mc_rtc::gui::NumberInput(
          "comZGainP", [this]() { return config().comZGainP; }, [this](double v) { config().comZGainP = v; }),
      mc_rtc::gui::NumberInput(
          "comZGainD", [this]() { return config().comZGainD; }, [this](double v) { config().comZGainD = v; }),
      mc_rtc::gui::Checkbox(
          "useTargetPoseForControlRobotAnchorFrame",
          [this]() { return config().useTargetPoseForControlRobotAnchorFrame; },
          [this]() {
            config().useTargetPoseForControlRobotAnchorFrame = !config().useTargetPoseForControlRobotAnchorFrame;
          }),
      mc_rtc::gui::Checkbox(
          "useActualComForWrenchDist", [this]() { return config().useActualComForWrenchDist; },
          [this]() { config().useActualComForWrenchDist = !config().useActualComForWrenchDist; }),
      mc_rtc::gui::ArrayInput(
          "actualComOffset", {"x", "y", "z"}, [this]() -> const Eigen::Vector3d & { return config().actualComOffset; },
          [this](const Eigen::Vector3d & v) { config().actualComOffset = v; }));

  gui.addElement({ctl().name(), config().name, "DcmEstimator"},
                 mc_rtc::gui::Checkbox(
                     "enableDcmEstimator", [this]() { return config().dcmEstimatorConfig.enableDcmEstimator; },
                     [this]() {
                       config().dcmEstimatorConfig.enableDcmEstimator = !config().dcmEstimatorConfig.enableDcmEstimator;
                     }),
                 mc_rtc::gui::ComboInput(
                     "dcmCorrectionMode", {"Bias", "Filter", "None"},
                     [this]() -> const std::string & { return config().dcmEstimatorConfig.dcmCorrectionMode; },
                     [this](const std::string & v) { config().dcmEstimatorConfig.dcmCorrectionMode = v; }),
                 mc_rtc::gui::NumberInput(
                     "biasDriftPerSecondStd", [this]() { return config().dcmEstimatorConfig.biasDriftPerSecondStd; },
                     [this](double v) {
                       config().dcmEstimatorConfig.biasDriftPerSecondStd = v;
                       dcmEstimator_->setBiasDriftPerSecond(v);
                     }),
                 mc_rtc::gui::NumberInput(
                     "dcmMeasureErrorStd", [this]() { return config().dcmEstimatorConfig.dcmMeasureErrorStd; },
                     [this](double v) {
                       config().dcmEstimatorConfig.dcmMeasureErrorStd = v;
                       dcmEstimator_->setDcmMeasureErrorStd(v);
                     }),
                 mc_rtc::gui::NumberInput(
                     "zmpMeasureErrorStd", [this]() { return config().dcmEstimatorConfig.zmpMeasureErrorStd; },
                     [this](double v) {
                       config().dcmEstimatorConfig.zmpMeasureErrorStd = v;
                       dcmEstimator_->setZmpMeasureErrorStd(v);
                     }),
                 mc_rtc::gui::ArrayInput(
                     "biasLimit", {"local X", "local Y"},
                     [this]() -> const Eigen::Vector2d & { return config().dcmEstimatorConfig.biasLimit; },
                     [this](const Eigen::Vector2d & v) {
                       config().dcmEstimatorConfig.biasLimit = v;
                       dcmEstimator_->setBiasLimit(v);
                     }),
                 mc_rtc::gui::ArrayLabel("estimatedLocalBias", {"local X", "local Y"},
                                         [this]() { return dcmEstimator_->getLocalBias(); }));

  gui.addElement(
      {ctl().name(), config().name, "Plot"}, mc_rtc::gui::ElementsStacking::Horizontal,
      mc_rtc::gui::Button(
          "Plot CoM-ZMP-X",
          [this, &gui]() {
            using namespace mc_rtc::gui;
            gui.addPlot(
                "CoM-ZMP-X", plot::X("t", [this]() { return ctl().t(); }),
                plot::Y(
                    "CoM_planned", [this]() { return ctl().comTask_->com().x(); }, Color::Blue, plot::Style::Dotted),
                plot::Y(
                    "CoM_controlRobot", [this]() { return ctl().robot().com().x(); }, Color::Green,
                    plot::Style::Dotted),
                plot::Y(
                    "CoM_realRobot", [this]() { return actualCom().x(); }, Color::Red, plot::Style::Dotted),
                plot::Y(
                    "CoM_realRobotUnbiased", [this]() { return actualComUnbiased().x(); },
                    mc_rtc::gui::Color(1.0, 0.5, 0.0), plot::Style::Dotted),
                plot::Y(
                    "ZMP_ref", [this]() { return refZmp_.x(); }, Color::Blue),
                plot::Y(
                    "ZMP_planned", [this]() { return plannedZmp_.x(); }, Color::Green),
                plot::Y(
                    "ZMP_control", [this]() { return controlZmp_.x(); }, Color::Magenta),
                plot::Y(
                    "ZMP_measured", [this]() { return measuredZMP_.x(); }, Color::Red),
                plot::Y(
                    "SupportRegion_min", [this]() { return supportRegion_[0].x(); }, Color::Black),
                plot::Y(
                    "SupportRegion_max", [this]() { return supportRegion_[1].x(); }, Color::Black));
          }),
      mc_rtc::gui::Button("Stop CoM-ZMP-X", [&gui]() { gui.removePlot("CoM-ZMP-X"); }));
  gui.addElement(
      {ctl().name(), config().name, "Plot"}, mc_rtc::gui::ElementsStacking::Horizontal,
      mc_rtc::gui::Button(
          "Plot CoM-ZMP-Y",
          [this, &gui]() {
            using namespace mc_rtc::gui;
            gui.addPlot(
                "CoM-ZMP-Y", plot::X("t", [this]() { return ctl().t(); }),
                plot::Y(
                    "CoM_planned", [this]() { return ctl().comTask_->com().y(); }, Color::Blue, plot::Style::Dotted),
                plot::Y(
                    "CoM_controlRobot", [this]() { return ctl().robot().com().y(); }, Color::Green,
                    plot::Style::Dotted),
                plot::Y(
                    "CoM_realRobot", [this]() { return actualCom().y(); }, Color::Red, plot::Style::Dotted),
                plot::Y(
                    "CoM_realRobotUnbiased", [this]() { return actualComUnbiased().y(); },
                    mc_rtc::gui::Color(1.0, 0.5, 0.0), plot::Style::Dotted),
                plot::Y(
                    "ZMP_ref", [this]() { return refZmp_.y(); }, Color::Blue),
                plot::Y(
                    "ZMP_planned", [this]() { return plannedZmp_.y(); }, Color::Green),
                plot::Y(
                    "ZMP_control", [this]() { return controlZmp_.y(); }, Color::Magenta),
                plot::Y(
                    "ZMP_measured", [this]() { return measuredZMP_.y(); }, Color::Red),
                plot::Y(
                    "SupportRegion_min", [this]() { return supportRegion_[0].y(); }, Color::Black),
                plot::Y(
                    "SupportRegion_max", [this]() { return supportRegion_[1].y(); }, Color::Black));
          }),
      mc_rtc::gui::Button("Stop CoM-ZMP-Y", [&gui]() { gui.removePlot("CoM-ZMP-Y"); }));
}

void CentroidalManager::removeFromGUI(mc_rtc::gui::StateBuilder & gui)
{
  gui.removeCategory({ctl().name(), config().name});
}

void CentroidalManager::addToLogger(mc_rtc::Logger & logger)
{
  logger.addLogEntry(config().name + "_Config_method", this, [this]() { return config().method; });
  logger.addLogEntry(config().name + "_Config_useActualStateForMpc", this,
                     [this]() { return config().useActualStateForMpc; });
  logger.addLogEntry(config().name + "_Config_enableZmpFeedback", this,
                     [this]() { return config().enableZmpFeedback; });
  logger.addLogEntry(config().name + "_Config_enableComZFeedback", this,
                     [this]() { return config().enableComZFeedback; });
  logger.addLogEntry(config().name + "_Config_dcmGainP", this, [this]() { return config().dcmGainP; });
  logger.addLogEntry(config().name + "_Config_zmpVelGain", this, [this]() { return config().zmpVelGain; });
  logger.addLogEntry(config().name + "_Config_comZGainP", this, [this]() { return config().comZGainP; });
  logger.addLogEntry(config().name + "_Config_comZGainD", this, [this]() { return config().comZGainD; });
  logger.addLogEntry(config().name + "_Config_refComZ", this, [this]() { return config().refComZ; });
  logger.addLogEntry(config().name + "_Config_useTargetPoseForControlRobotAnchorFrame", this,
                     [this]() { return config().useTargetPoseForControlRobotAnchorFrame; });
  logger.addLogEntry(config().name + "_Config_useActualComForWrenchDist", this,
                     [this]() { return config().useActualComForWrenchDist; });
  logger.addLogEntry(config().name + "_Config_actualComOffset", this, [this]() { return config().actualComOffset; });

  MC_RTC_LOG_HELPER(config().name + "_CoM_MPC", mpcCom_);
  logger.addLogEntry(config().name + "_CoM_planned", this, [this]() { return ctl().comTask_->com(); });
  logger.addLogEntry(config().name + "_CoM_controlRobot", this, [this]() { return ctl().robot().com(); });
  logger.addLogEntry(config().name + "_CoM_realRobot", this, [this]() { return actualCom(); });
  logger.addLogEntry(config().name + "_CoM_realRobotUnbiased", this, [this]() { return actualComUnbiased(); });

  MC_RTC_LOG_HELPER(config().name + "_forceZ_planned", plannedForceZ_);
  MC_RTC_LOG_HELPER(config().name + "_forceZ_control", controlForceZ_);

  logger.addLogEntry(config().name + "_ZMP_ref", this, [this]() { return refZmp_; });
  MC_RTC_LOG_HELPER(config().name + "_ZMP_planned", plannedZmp_);
  MC_RTC_LOG_HELPER(config().name + "_ZMP_control", controlZmp_);
  logger.addLogEntry(config().name + "_ZMP_controlWrenchDist", this, [this]() {
    return wrenchDist_ ? calcZmp(ForceColl::calcWrenchList(contactList_, wrenchDist_->resultWrenchRatio_), refZmp_.z())
                       : Eigen::Vector3d::Zero();
  });
  MC_RTC_LOG_HELPER(config().name + "_ZMP_measured", measuredZMP_);
  logger.addLogEntry(config().name + "_ZMP_SupportRegion_min", this, [this]() { return supportRegion_[0]; });
  logger.addLogEntry(config().name + "_ZMP_SupportRegion_max", this, [this]() { return supportRegion_[1]; });

  logger.addLogEntry(config().name + "_CentroidalMomentum_controlRobot", this, [this]() {
    return rbd::computeCentroidalMomentum(ctl().robot().mb(), ctl().robot().mbc(), ctl().robot().com());
  });

  logger.addLogEntry(config().name + "_DcmEstimator_enableDcmEstimator", this,
                     [this]() { return config().dcmEstimatorConfig.enableDcmEstimator; });
  logger.addLogEntry(config().name + "_DcmEstimator_dcmCorrectionMode", this,
                     [this]() { return config().dcmEstimatorConfig.dcmCorrectionMode; });
  logger.addLogEntry(config().name + "_DcmEstimator_localBias", this,
                     [this]() { return dcmEstimator_->getLocalBias(); });
}

void CentroidalManager::removeFromLogger(mc_rtc::Logger & logger)
{
  logger.removeLogEntries(this);
}

bool CentroidalManager::setRefComZ(double refComZ, double startTime, double interpDuration)
{
  if(startTime < ctl().t())
  {
    mc_rtc::log::warning("[CentroidalManager] Ignore reference CoM Z position with past time: {} < {}", startTime,
                         ctl().t());
    return false;
  }
  auto it = std::next(refComZFunc_->points().rbegin());
  if(startTime < it->first)
  {
    mc_rtc::log::warning("[CentroidalManager] Ignore reference CoM Z position with time before the existing "
                         "interpolation points: {} < {}",
                         startTime, it->first);
    return false;
  }

  refComZFunc_->appendPoint(std::make_pair(startTime, it->second));
  refComZFunc_->appendPoint(std::make_pair(startTime + interpDuration, refComZ));
  refComZFunc_->calcCoeff();

  return true;
}

void CentroidalManager::setAnchorFrame()
{
  std::string anchorName = "KinematicAnchorFrame::" + ctl().robot().name();
  if(ctl().datastore().has(anchorName))
  {
    ctl().datastore().remove(anchorName);
  }
  ctl().datastore().make_call(anchorName, [this](const mc_rbdyn::Robot & robot) { return calcAnchorFrame(robot); });
}

double CentroidalManager::calcRefComZ(double t, int derivOrder) const
{
  if(derivOrder == 0)
  {
    return (*refComZFunc_)(t);
  }
  else
  {
    return refComZFunc_->derivative(t, derivOrder);
  }
}

double CentroidalManager::getfootSurfaceDiff() const
{
  double footSurfaceDiff = 0.0;
  int contactCount = 0;

  // Get contact feet
  contactCount = static_cast<int>(ctl().footManager_->getCurrentContactFeet().size());

  // Is Robot on two feet?
  if(contactCount == 2 || isDangerFloor == true)
  {
    const auto & leftSurfacePose = ctl().robot().surfacePose(ctl().footManager_->surfaceName(Foot::Left));
    const auto & rightSurfacePose = ctl().robot().surfacePose(ctl().footManager_->surfaceName(Foot::Right));
    footSurfaceDiff = std::abs(leftSurfacePose.translation().z() - rightSurfacePose.translation().z());
    return footSurfaceDiff;
  }

  return 0.0;
}

double CentroidalManager::compliantFloorCorrection(double footSurfaceDiff) const
{
  double comZ = actualCom().z();                         // 実際のCoM高さ
  double comVelZ = ctl().realRobot().comVelocity().z();  // CoM速度Z

  double plannedGroundZ = ctl().footManager_->calcRefGroundPosZ(ctl().t());
  double supportZ = plannedGroundZ - footSurfaceDiff;
  double relZ = comZ - supportZ;      // 床からの相対高さ
  double relVelZ = comVelZ;

  // 床の物理定数
  const double K = config().floorSpringK;
  const double D = config().floorDampingD;

  // a = F/m = (K*Δz + D*Δẋ)/m
  return (K * (-relZ) + D * (-relVelZ)) / robotMass_;
}

double CentroidalManager::estimateSupportSink() const
{
  // 実足の平均支持高さ
  double supportZ = 0.0;
  int n = 0;
  for(const auto & foot : ctl().footManager_->getCurrentContactFeet())
  {
    supportZ += ctl().robot().surfacePose(ctl().footManager_->surfaceName(foot)).translation().z();
    ++n;
  }
  if(n == 0) { return 0.0; }
  supportZ /= static_cast<double>(n);

  // 計画上の地面高さ（剛床前提の参照）
  double plannedGroundZ = ctl().footManager_->calcRefGroundPosZ(ctl().t());
  config().lastSupportZ_ = supportZ;

  // 正なら「計画に対して実際の支持面が沈んだ」量
  return plannedGroundZ - supportZ;
}

sva::PTransformd CentroidalManager::calcAnchorFrame(const mc_rbdyn::Robot & robot) const
{
  double leftFootSupportRatio = ctl().footManager_->leftFootSupportRatio();
  bool isControlRobot = (&(ctl().robot()) == &robot);

  if(isControlRobot && config().useTargetPoseForControlRobotAnchorFrame)
  {
    return sva::interpolate(ctl().footManager_->targetFootPose(Foot::Right),
                            ctl().footManager_->targetFootPose(Foot::Left), leftFootSupportRatio);
  }
  else
  {
    return sva::interpolate(robot.surfacePose(ctl().footManager_->surfaceName(Foot::Right)),
                            robot.surfacePose(ctl().footManager_->surfaceName(Foot::Left)), leftFootSupportRatio);
  }
}

Eigen::Vector3d CentroidalManager::actualCom() const
{
  return ctl().realRobot().com() + config().actualComOffset;
}

Eigen::Vector3d CentroidalManager::actualComUnbiased() const
{
  Eigen::Vector3d com = actualCom();
  com.head<2>() -= dcmEstimator_->getBias();
  return com;
}

Eigen::Vector3d CentroidalManager::calcZmp(const std::unordered_map<Foot, sva::ForceVecd> & wrenchList,
                                           double zmpPlaneHeight,
                                           const Eigen::Vector3d & zmpPlaneNormal) const
{
  sva::ForceVecd totalWrench = sva::ForceVecd::Zero();

  // モーメントの総和
  for(const auto & wrenchKV : wrenchList)
  {
    totalWrench += wrenchKV.second;
  }

  Eigen::Vector3d zmpPlaneOrigin = Eigen::Vector3d(0, 0, zmpPlaneHeight);
  Eigen::Vector3d zmp = zmpPlaneOrigin;

  if(totalWrench.force().z() > 0)
  {
    Eigen::Vector3d momentInZmpPlane = totalWrench.moment() - zmpPlaneOrigin.cross(totalWrench.force());
    zmp += zmpPlaneNormal.cross(momentInZmpPlane) / totalWrench.force().z();
  }

  return zmp;
}

Eigen::Vector3d CentroidalManager::calcPlannedComAccel() const
{
  Eigen::Vector3d plannedComAccel;
  plannedComAccel << plannedForceZ_ / (robotMass_ * (mpcCom_.z() - refZmp_.z()))
                         * (mpcCom_.head<2>() - plannedZmp_.head<2>()),
      plannedForceZ_ / robotMass_;
  plannedComAccel.z() -= CCC::constants::g;
  return plannedComAccel;
}
