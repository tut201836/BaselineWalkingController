#include <algorithm>
#include <limits>

#include <mc_filter/utils/clamp.h>
#include <mc_rtc/gui/ArrayInput.h>
#include <mc_rtc/gui/Checkbox.h>
#include <mc_rtc/gui/ComboInput.h>
#include <mc_rtc/gui/IntegerInput.h>
#include <mc_rtc/gui/Label.h>
#include <mc_rtc/gui/NumberInput.h>
#include <mc_rtc/gui/Polygon.h>
#include <mc_tasks/FirstOrderImpedanceTask.h>
#include <mc_tasks/OrientationTask.h>

#include <ForceColl/Contact.h>

#include <BaselineWalkingController/BaselineWalkingController.h>
#include <BaselineWalkingController/FootManager.h>
#include <BaselineWalkingController/MathUtils.h>
#include <BaselineWalkingController/swing/SwingTrajCubicSplineSimple.h>
#include <BaselineWalkingController/swing/SwingTrajIndHorizontalVertical.h>
#include <BaselineWalkingController/swing/SwingTrajLandingSearch.h>
#include <BaselineWalkingController/swing/SwingTrajVariableTaskGain.h>

using namespace BWC;

void FootManager::Configuration::load(const mc_rtc::Configuration & mcRtcConfig)
{
  mcRtcConfig("name", name);
  mcRtcConfig("footstepDuration", footstepDuration);
  mcRtcConfig("doubleSupportRatio", doubleSupportRatio);
  if(mcRtcConfig.has("deltaTransLimit"))
  {
    deltaTransLimit = mcRtcConfig("deltaTransLimit");
    deltaTransLimit[2] = mc_rtc::constants::toRad(deltaTransLimit[2]);
  }
  if(mcRtcConfig.has("midToFootTranss"))
  {
    for(const auto & foot : Feet::Both)
    {
      mcRtcConfig("midToFootTranss")(std::to_string(foot), midToFootTranss.at(foot));
    }
  }
  if(mcRtcConfig.has("footTaskGain"))
  {
    footTaskGain = TaskGain(mcRtcConfig("footTaskGain"));
  }
  mcRtcConfig("zmpHorizon", zmpHorizon);
  mcRtcConfig("zmpOffset", zmpOffset);
  mcRtcConfig("defaultSwingTrajType", defaultSwingTrajType);
  mcRtcConfig("overwriteLandingPose", overwriteLandingPose);
  mcRtcConfig("stopSwingTrajForTouchDownFoot", stopSwingTrajForTouchDownFoot);
  mcRtcConfig("keepPoseForTouchDownFoot", keepPoseForTouchDownFoot);
  mcRtcConfig("enableWrenchDistForTouchDownFoot", enableWrenchDistForTouchDownFoot);
  mcRtcConfig("enableArmSwing", enableArmSwing);
  mcRtcConfig("fricCoeff", fricCoeff);
  mcRtcConfig("touchDownRemainingDuration", touchDownRemainingDuration);
  mcRtcConfig("touchDownPosError", touchDownPosError);
  mcRtcConfig("touchDownForceZ", touchDownForceZ);
  if(mcRtcConfig.has("impedanceGains"))
  {
    mcRtcConfig("impedanceGains")("SingleSupport", impGains.at("SingleSupport"));
    mcRtcConfig("impedanceGains")("DoubleSupport", impGains.at("DoubleSupport"));
    mcRtcConfig("impedanceGains")("Swing", impGains.at("Swing"));
  }
  if(mcRtcConfig.has("jointAnglesForArmSwing"))
  {
    mcRtcConfig("jointAnglesForArmSwing")("Nominal", jointAnglesForArmSwing.at("Nominal"));
    mcRtcConfig("jointAnglesForArmSwing")("Left", jointAnglesForArmSwing.at("Left"));
    mcRtcConfig("jointAnglesForArmSwing")("Right", jointAnglesForArmSwing.at("Right"));
  }
}

void FootManager::VelModeData::Configuration::load(const mc_rtc::Configuration & mcRtcConfig)
{
  if(mcRtcConfig.has("footstepQueueSize"))
  {
    footstepQueueSize = mcRtcConfig("footstepQueueSize");
    if(footstepQueueSize < 3)
    {
      footstepQueueSize = 3;
      mc_rtc::log::warning("[FootManager::VelModeData] footstepQueueSize must be at least 3.");
    }
  }
  mcRtcConfig("enableOnlineFootstepUpdate", enableOnlineFootstepUpdate);
}

void FootManager::VelModeData::reset(bool enabled)
{
  enabled_ = enabled;
  targetVel_.setZero();
}

FootManager::FootManager(BaselineWalkingController * ctlPtr, const mc_rtc::Configuration & mcRtcConfig)
: ctlPtr_(ctlPtr), zmpFunc_(std::make_shared<TrajColl::CubicInterpolator<Eigen::Vector3d>>()),
  groundPosZFunc_(std::make_shared<TrajColl::CubicInterpolator<double>>()),
  baseYawFunc_(std::make_shared<TrajColl::CubicInterpolator<Eigen::Matrix3d, Eigen::Vector3d>>())
{
  config_.load(mcRtcConfig);

  if(mcRtcConfig.has("VelMode"))
  {
    velModeData_.config_.load(mcRtcConfig("VelMode"));
  }

  if(mcRtcConfig.has("SwingTraj"))
  {
    SwingTrajCubicSplineSimple::loadDefaultConfig(
        mcRtcConfig("SwingTraj")("CubicSplineSimple", mc_rtc::Configuration{}));
    SwingTrajIndHorizontalVertical::loadDefaultConfig(
        mcRtcConfig("SwingTraj")("IndHorizontalVertical", mc_rtc::Configuration{}));
    SwingTrajVariableTaskGain::loadDefaultConfig(mcRtcConfig("SwingTraj")("VariableTaskGain", mc_rtc::Configuration{}));
    SwingTrajLandingSearch::loadDefaultConfig(mcRtcConfig("SwingTraj")("LandingSearch", mc_rtc::Configuration{}));
  }
}

void FootManager::reset()
{
  velModeData_.reset(false);

  footstepQueue_.clear();
  prevFootstep_.reset();

  targetFootPoses_.clear();
  targetFootVels_.clear();
  targetFootAccels_.clear();
  footTaskGains_.clear();
  trajStartFootPoseFuncs_.clear();
  for(const auto & foot : Feet::Both)
  {
    targetFootPoses_.emplace(foot, ctl().robot().surfacePose(surfaceName(foot)));
    targetFootVels_.emplace(foot, sva::MotionVecd::Zero());
    targetFootAccels_.emplace(foot, sva::MotionVecd::Zero());
    footTaskGains_.emplace(foot, config_.footTaskGain);
    trajStartFootPoseFuncs_.emplace(foot, nullptr);
  }
  trajStartFootPoses_ = targetFootPoses_;

  supportPhase_ = SupportPhase::DoubleSupport;

  Eigen::Vector3d targetZmp = calcZmpWithOffset(targetFootPoses_);
  zmpFunc_->clearPoints();
  zmpFunc_->appendPoint(std::make_pair(ctl().t(), targetZmp));
  zmpFunc_->appendPoint(std::make_pair(ctl().t() + config_.zmpHorizon, targetZmp));
  zmpFunc_->calcCoeff();

  double refGroundPosZ =
      0.5 * (targetFootPoses_.at(Foot::Left).translation().z() + targetFootPoses_.at(Foot::Right).translation().z());
  groundPosZFunc_->clearPoints();
  groundPosZFunc_->appendPoint(std::make_pair(ctl().t(), refGroundPosZ));
  groundPosZFunc_->appendPoint(std::make_pair(ctl().t() + config_.zmpHorizon, refGroundPosZ));
  groundPosZFunc_->calcCoeff();

  contactFootPosesList_.emplace(ctl().t(), targetFootPoses_);

  swingFootstep_ = nullptr;
  swingTraj_.reset();

  baseYawFunc_->clearPoints();

  armSwingFunc_.reset();

  touchDown_ = false;

  for(const auto & foot : Feet::Both)
  {
    impGainTypes_.emplace(foot, "DoubleSupport");

    // init landing Z offsets
    landingZOffset_.at(foot) = 0.0;
  }

  requireImpGainUpdate_ = true;

  if(config_.jointAnglesForArmSwing.at("Nominal").empty() && config_.jointAnglesForArmSwing.at("Left").size() > 0)
  {
    auto postureTask = ctl().getPostureTask(ctl().robot().name());
    for(const auto & jointAngleKV : config_.jointAnglesForArmSwing.at("Left"))
    {
      config_.jointAnglesForArmSwing.at("Nominal")[jointAngleKV.first] =
          postureTask->posture()[ctl().robot().jointIndexByName(jointAngleKV.first)];
    }
  }

  // 追加：世界基準Zを記録
  worldGroundZ0_ = refGroundPosZ;

  groundPosZFunc_->clearPoints();
  groundPosZFunc_->appendPoint(std::make_pair(ctl().t(), refGroundPosZ));
  groundPosZFunc_->appendPoint(std::make_pair(ctl().t() + config_.zmpHorizon, refGroundPosZ));
  groundPosZFunc_->calcCoeff();
}

void FootManager::update()
{
  updateFootTraj();
  updateZmpTraj();
  if(velModeData_.enabled_)
  {
    updateVelMode();
  }
}

void FootManager::stop()
{
  removeFromGUI(*ctl().gui());
  removeFromLogger(ctl().logger());
}

void FootManager::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  gui.addElement({ctl().name(), config_.name, "Status"},
                 mc_rtc::gui::Label("supportPhase", [this]() { return std::to_string(supportPhase_); }),
                 mc_rtc::gui::Label("footstepQueueSize", [this]() { return std::to_string(footstepQueue_.size()); }));
  gui.addElement({ctl().name(), config_.name, "Status"}, mc_rtc::gui::ElementsStacking::Horizontal,
                 mc_rtc::gui::Label("LeftFootSurface", [this]() { return surfaceName(Foot::Left); }),
                 mc_rtc::gui::Label("RightFootSurface", [this]() { return surfaceName(Foot::Right); }));
  gui.addElement({ctl().name(), config_.name, "Status"}, mc_rtc::gui::ElementsStacking::Horizontal,
                 mc_rtc::gui::Label("LeftImpGainType", [this]() { return impGainTypes_.at(Foot::Left); }),
                 mc_rtc::gui::Label("RightImpGainType", [this]() { return impGainTypes_.at(Foot::Right); }));

  gui.addElement(
      {ctl().name(), config_.name, "Config"},
      mc_rtc::gui::NumberInput(
          "footstepDuration", [this]() { return config_.footstepDuration; },
          [this](double v) { config_.footstepDuration = v; }),
      mc_rtc::gui::NumberInput(
          "doubleSupportRatio", [this]() { return config_.doubleSupportRatio; },
          [this](double v) { config_.doubleSupportRatio = v; }),
      mc_rtc::gui::ArrayInput(
          "deltaTransLimit", {"x", "y", "theta"},
          [this]() -> Eigen::Vector3d {
            return Eigen::Vector3d(config_.deltaTransLimit[0], config_.deltaTransLimit[1],
                                   mc_rtc::constants::toDeg(config_.deltaTransLimit[2]));
          },
          [this](const Eigen::Vector3d & v) {
            config_.deltaTransLimit = Eigen::Vector3d(v[0], v[1], mc_rtc::constants::toRad(v[2]));
          }),
      mc_rtc::gui::ArrayInput(
          "footTaskStiffness", {"ax", "ay", "az", "tx", "ty", "tz"},
          [this]() -> const sva::MotionVecd & { return config_.footTaskGain.stiffness; },
          [this](const Eigen::Vector6d & v) { config_.footTaskGain.stiffness = sva::MotionVecd(v); }),
      mc_rtc::gui::ArrayInput(
          "footTaskDamping", {"ax", "ay", "az", "tx", "ty", "tz"},
          [this]() -> const sva::MotionVecd & { return config_.footTaskGain.damping; },
          [this](const Eigen::Vector6d & v) { config_.footTaskGain.damping = sva::MotionVecd(v); }),
      mc_rtc::gui::ArrayInput(
          "zmpOffset", {"x", "y", "z"}, [this]() -> const Eigen::Vector3d & { return config_.zmpOffset; },
          [this](const Eigen::Vector3d & v) { config_.zmpOffset = v; }),
      mc_rtc::gui::ComboInput(
          "defaultSwingTrajType", {"CubicSplineSimple", "IndHorizontalVertical", "VariableTaskGain", "LandingSearch"},
          [this]() { return config_.defaultSwingTrajType; },
          [this](const std::string & v) { config_.defaultSwingTrajType = v; }),
      mc_rtc::gui::Checkbox(
          "overwriteLandingPose", [this]() { return config_.overwriteLandingPose; },
          [this]() { config_.overwriteLandingPose = !config_.overwriteLandingPose; }),
      mc_rtc::gui::Checkbox(
          "stopSwingTrajForTouchDownFoot", [this]() { return config_.stopSwingTrajForTouchDownFoot; },
          [this]() { config_.stopSwingTrajForTouchDownFoot = !config_.stopSwingTrajForTouchDownFoot; }),
      mc_rtc::gui::Checkbox(
          "keepPoseForTouchDownFoot", [this]() { return config_.keepPoseForTouchDownFoot; },
          [this]() { config_.keepPoseForTouchDownFoot = !config_.keepPoseForTouchDownFoot; }),
      mc_rtc::gui::Checkbox(
          "enableWrenchDistForTouchDownFoot", [this]() { return config_.enableWrenchDistForTouchDownFoot; },
          [this]() { config_.enableWrenchDistForTouchDownFoot = !config_.enableWrenchDistForTouchDownFoot; }),
      mc_rtc::gui::Checkbox(
          "enableArmSwing", [this]() { return config_.enableArmSwing; },
          [this]() { config_.enableArmSwing = !config_.enableArmSwing; }),
      mc_rtc::gui::NumberInput(
          "fricCoeff", [this]() { return config_.fricCoeff; }, [this](double v) { config_.fricCoeff = v; }),
      mc_rtc::gui::NumberInput(
          "touchDownRemainingDuration", [this]() { return config_.touchDownRemainingDuration; },
          [this](double v) { config_.touchDownRemainingDuration = v; }),
      mc_rtc::gui::NumberInput(
          "touchDownPosError", [this]() { return config_.touchDownPosError; },
          [this](double v) { config_.touchDownPosError = v; }),
      mc_rtc::gui::NumberInput(
          "touchDownForceZ", [this]() { return config_.touchDownForceZ; },
          [this](double v) { config_.touchDownForceZ = v; }),
      mc_rtc::gui::Label("jointsForArmSwing", [this]() {
        std::string s;
        for(const auto & jointAngleKV : config_.jointAnglesForArmSwing.at("Nominal"))
        {
          if(!s.empty())
          {
            s += " / ";
          }
          s += jointAngleKV.first;
        }
        return s;
      }));

  gui.addElement(
      {ctl().name(), config_.name, "Config", "VelMode"},
      mc_rtc::gui::IntegerInput(
          "footstepQueueSize", [this]() { return velModeData_.config_.footstepQueueSize; },
          [this](int footstepQueueSize) { velModeData_.config_.footstepQueueSize = std::max(footstepQueueSize, 3); }),
      mc_rtc::gui::Checkbox(
          "enableOnlineFootstepUpdate", [this]() { return velModeData_.config_.enableOnlineFootstepUpdate; },
          [this]() {
            velModeData_.config_.enableOnlineFootstepUpdate = !velModeData_.config_.enableOnlineFootstepUpdate;
          }));

  for(const auto & impGainKV : config_.impGains)
  {
    const auto & impGainType = impGainKV.first;
    gui.addElement({ctl().name(), config_.name, "ImpedanceGains", impGainType},
                   mc_rtc::gui::ArrayInput(
                       "Damper", {"cx", "cy", "cz", "fx", "fy", "fz"},
                       [this, impGainType]() -> const sva::ImpedanceVecd & {
                         return config_.impGains.at(impGainType).damper().vec();
                       },
                       [this, impGainType](const Eigen::Vector6d & v) {
                         config_.impGains.at(impGainType).damper().vec(v);
                         requireImpGainUpdate_ = true;
                       }));
    gui.addElement({ctl().name(), config_.name, "ImpedanceGains", impGainType},
                   mc_rtc::gui::ArrayInput(
                       "Spring", {"cx", "cy", "cz", "fx", "fy", "fz"},
                       [this, impGainType]() -> const sva::ImpedanceVecd & {
                         return config_.impGains.at(impGainType).spring().vec();
                       },
                       [this, impGainType](const Eigen::Vector6d & v) {
                         config_.impGains.at(impGainType).spring().vec(v);
                         requireImpGainUpdate_ = true;
                       }));
    gui.addElement({ctl().name(), config_.name, "ImpedanceGains", impGainType},
                   mc_rtc::gui::ArrayInput(
                       "Wrench", {"cx", "cy", "cz", "fx", "fy", "fz"},
                       [this, impGainType]() -> const sva::ImpedanceVecd & {
                         return config_.impGains.at(impGainType).wrench().vec();
                       },
                       [this, impGainType](const Eigen::Vector6d & v) {
                         config_.impGains.at(impGainType).wrench().vec(v);
                         requireImpGainUpdate_ = true;
                       }));
  }

  SwingTrajCubicSplineSimple::addConfigToGUI(gui, {ctl().name(), "SwingTraj", "CubicSplineSimple"});
  SwingTrajIndHorizontalVertical::addConfigToGUI(gui, {ctl().name(), "SwingTraj", "IndHorizontalVertical"});
  SwingTrajVariableTaskGain::addConfigToGUI(gui, {ctl().name(), "SwingTraj", "VariableTaskGain"});
  SwingTrajLandingSearch::addConfigToGUI(gui, {ctl().name(), "SwingTraj", "LandingSearch"});
}

void FootManager::removeFromGUI(mc_rtc::gui::StateBuilder & gui)
{
  gui.removeCategory({ctl().name(), config_.name});

  SwingTrajCubicSplineSimple::removeConfigFromGUI(gui, {ctl().name(), "SwingTraj", "CubicSplineSimple"});
  SwingTrajIndHorizontalVertical::removeConfigFromGUI(gui, {ctl().name(), "SwingTraj", "IndHorizontalVertical"});
  SwingTrajVariableTaskGain::removeConfigFromGUI(gui, {ctl().name(), "SwingTraj", "VariableTaskGain"});
  SwingTrajLandingSearch::removeConfigFromGUI(gui, {ctl().name(), "SwingTraj", "LandingSearch"});
}

void FootManager::addToLogger(mc_rtc::Logger & logger)
{
  logger.addLogEntry(config_.name + "_footstepQueueSize", this, [this]() { return footstepQueue_.size(); });
  logger.addLogEntry(config_.name + "_walking", this, [this]() {
    return static_cast<int>(!footstepQueue_.empty() && footstepQueue_.front().transitStartTime <= ctl().t());
  });

  for(const auto & foot : Feet::Both)
  {
    logger.addLogEntry(config_.name + "_targetFootPose_" + std::to_string(foot), this,
                       [this, foot]() -> const sva::PTransformd & { return targetFootPoses_.at(foot); });

    logger.addLogEntry(config_.name + "_targetFootVel_" + std::to_string(foot), this,
                       [this, foot]() -> const sva::MotionVecd & { return targetFootVels_.at(foot); });

    logger.addLogEntry(config_.name + "_targetFootAccel_" + std::to_string(foot), this,
                       [this, foot]() -> const sva::MotionVecd & { return targetFootAccels_.at(foot); });

    logger.addLogEntry(config_.name + "_footTaskStiffness_" + std::to_string(foot), this,
                       [this, foot]() -> const sva::MotionVecd & { return footTaskGains_.at(foot).stiffness; });

    logger.addLogEntry(config_.name + "_footTaskDamping_" + std::to_string(foot), this,
                       [this, foot]() -> const sva::MotionVecd & { return footTaskGains_.at(foot).damping; });
  }
  logger.addLogEntry(config_.name + "_swingTrajType", this,
                     [this]() -> std::string { return swingTraj_ ? swingTraj_->type() : "None"; });

  logger.addLogEntry(config_.name + "_supportPhase", this, [this]() { return std::to_string(supportPhase_); });

  logger.addLogEntry(config_.name + "_refZmp", this, [this]() { return calcRefZmp(ctl().t()); });

  logger.addLogEntry(config_.name + "_refGroundPosZ", this, [this]() { return calcRefGroundPosZ(ctl().t()); });

  logger.addLogEntry(config_.name + "_leftFootSupportRatio", this, [this]() { return leftFootSupportRatio(); });

  logger.addLogEntry(config_.name + "_velMode", this,
                     [this]() -> std::string { return velModeData_.enabled_ ? "ON" : "OFF"; });
  logger.addLogEntry(config_.name + "_targetVel", this, [this]() { return velModeData_.targetVel_; });

  logger.addLogEntry(config_.name + "_touchDown", this, [this]() { return touchDown_; });

  logger.addLogEntry(config_.name + "_touchDownRemainingDuration", this,
                     [this]() { return touchDownRemainingDuration(); });

  for(const auto & foot : Feet::Both)
  {
    logger.addLogEntry(config_.name + "_impGainType_" + std::to_string(foot), this,
                       [this, foot]() -> const std::string & { return impGainTypes_.at(foot); });
  }
}

void FootManager::removeFromLogger(mc_rtc::Logger & logger)
{
  logger.removeLogEntries(this);
}

const std::string & FootManager::surfaceName(const Foot & foot) const
{
  return ctl().footTasks_.at(foot)->surface();
}

Footstep FootManager::makeFootstep(const Foot & foot,
                                   const sva::PTransformd & footMidpose,
                                   double startTime,
                                   const mc_rtc::Configuration & swingTrajConfig) const
{
  return Footstep(foot, config_.midToFootTranss.at(foot) * footMidpose, startTime,
                  startTime + 0.5 * config_.doubleSupportRatio * config_.footstepDuration,
                  startTime + (1.0 - 0.5 * config_.doubleSupportRatio) * config_.footstepDuration,
                  startTime + config_.footstepDuration, swingTrajConfig);
}

bool FootManager::appendFootstep(const Footstep & newFootstep)
{
  // Check time of new footstep
  if(newFootstep.transitStartTime < ctl().t())
  {
    mc_rtc::log::error("[FootManager] Ignore a new footstep with past time: {} < {}", newFootstep.transitStartTime,
                       ctl().t());
    return false;
  }
  if(!footstepQueue_.empty())
  {
    const Footstep & lastFootstep = footstepQueue_.back();
    if(newFootstep.transitStartTime < lastFootstep.transitEndTime)
    {
      mc_rtc::log::error("[FootManager] Ignore a new footstep earlier than the last footstep: {} < {}",
                         newFootstep.transitStartTime, lastFootstep.transitEndTime);
      return false;
    }
  }

  // Push to the queue
  footstepQueue_.push_back(newFootstep);

  return true;
}

void FootManager::clearFootstepQueue()
{
  if(footstepQueue_.empty())
  {
    return;
  }

  if(swingFootstep_ == &(footstepQueue_.front()))
  {
    if(footstepQueue_.size() >= 2)
    {
      footstepQueue_.erase(footstepQueue_.begin() + 1, footstepQueue_.end());
    }
  }
  else
  {
    footstepQueue_.clear();
  }
}

Eigen::Vector3d FootManager::clampDeltaTrans(const Eigen::Vector3d & deltaTrans, const Foot & foot)
{
  Eigen::Vector3d deltaTransMax = config_.deltaTransLimit;
  Eigen::Vector3d deltaTransMin = -1 * config_.deltaTransLimit;
  if(foot == Foot::Left)
  {
    deltaTransMin.y() = 0;
  }
  else
  {
    deltaTransMax.y() = 0;
  }
  return mc_filter::utils::clamp(deltaTrans, deltaTransMin, deltaTransMax);
}

Eigen::Vector3d FootManager::calcRefZmp(double t, int derivOrder) const
{
  if(derivOrder == 0)
  {
    return (*zmpFunc_)(t);
  }
  else
  {
    return zmpFunc_->derivative(t, derivOrder);
  }
}

double FootManager::calcRefGroundPosZ(double t, int derivOrder) const
{
  if(derivOrder == 0)
  {
    return (*groundPosZFunc_)(t);
  }
  else
  {
    return groundPosZFunc_->derivative(t, derivOrder);
  }
}

std::unordered_map<Foot, sva::PTransformd> FootManager::calcContactFootPoses(double t) const
{
  auto it = contactFootPosesList_.upper_bound(t);
  if(it == contactFootPosesList_.begin())
  {
    return std::unordered_map<Foot, sva::PTransformd>{};
  }
  else
  {
    it--;
    return it->second;
  }
}

std::set<Foot> FootManager::getCurrentContactFeet() const
{
  if(supportPhase_ == SupportPhase::DoubleSupport)
  {
    return Feet::Both;
  }
  else
  {
    if(config_.enableWrenchDistForTouchDownFoot && touchDown_)
    {
      return Feet::Both;
    }
    else
    {
      if(supportPhase_ == SupportPhase::LeftSupport)
      {
        return std::set<Foot>{Foot::Left};
      }
      else // if(supportPhase_ == SupportPhase::RightSupport)
      {
        return std::set<Foot>{Foot::Right};
      }
    }
  }
}

std::unordered_map<Foot, std::shared_ptr<ForceColl::Contact>> FootManager::calcCurrentContactList() const
{
  // Set contactList
  std::unordered_map<Foot, std::shared_ptr<ForceColl::Contact>> contactList;
  for(const auto & foot : getCurrentContactFeet())
  {
    const auto & surface = ctl().robot().surface(surfaceName(foot));
    contactList.emplace(
        foot, std::make_shared<ForceColl::SurfaceContact>(std::to_string(foot), config_.fricCoeff,
                                                          calcSurfaceVertexList(surface, sva::PTransformd::Identity()),
                                                          targetFootPoses_.at(foot)));
  }

  return contactList;
}

double FootManager::leftFootSupportRatio() const
{
  if(footstepQueue_.empty())
  {
    return 0.5;
  }

  const Footstep & footstep = footstepQueue_.front();

  if(ctl().t() <= footstep.transitStartTime)
  {
    return 0.5;
  }
  else if(ctl().t() <= footstep.swingStartTime)
  {
    double ratio = (ctl().t() - footstep.transitStartTime) / (footstep.swingStartTime - footstep.transitStartTime);
    return -0.5 * sign(footstep.foot) * mc_filter::utils::clamp(ratio, 0.0, 1.0) + 0.5;
  }
  else if(ctl().t() <= footstep.swingEndTime)
  {
    return footstep.foot == Foot::Left ? 0 : 1;
  }
  else // if(ctl().t() <= footstep.transitEndTime)
  {
    double ratio = (footstep.transitEndTime - ctl().t()) / (footstep.transitEndTime - footstep.swingEndTime);
    return -0.5 * sign(footstep.foot) * mc_filter::utils::clamp(ratio, 0.0, 1.0) + 0.5;
  }
}

Eigen::Vector3d FootManager::calcZmpWithOffset(const Foot & foot, const sva::PTransformd & footPose) const
{
  Eigen::Vector3d zmpOffset = config_.zmpOffset;
  if(foot == Foot::Right)
  {
    zmpOffset.y() *= -1;
  }
  return (sva::PTransformd(zmpOffset) * footPose).translation();
}

Eigen::Vector3d FootManager::calcZmpWithOffset(const std::unordered_map<Foot, sva::PTransformd> & footPoses) const
{
  if(footPoses.size() == 0)
  {
    mc_rtc::log::error("[FootManager] footPoses is empty in calcZmpWithOffset.");
    return Eigen::Vector3d::Zero();
  }
  else if(footPoses.size() == 1)
  {
    return calcZmpWithOffset(footPoses.begin()->first, footPoses.begin()->second);
  }
  else // if(footPoses.size() == 2)
  {
    return 0.5
           * (calcZmpWithOffset(Foot::Left, footPoses.at(Foot::Left))
              + calcZmpWithOffset(Foot::Right, footPoses.at(Foot::Right)));
  }
}

bool FootManager::walkToRelativePose(const Eigen::Vector3d & targetTrans,
                                     int lastFootstepNum,
                                     const std::vector<Eigen::Vector3d> & waypointTransList)
{
  if(footstepQueue_.size() > 0)
  {
    mc_rtc::log::error("[FootManager] walkToRelativePose is available only when the footstep queue is empty: {}",
                       footstepQueue_.size());
    return false;
  }

  auto convertTo2d = [](const sva::PTransformd & pose) -> Eigen::Vector3d {
    return Eigen::Vector3d(pose.translation().x(), pose.translation().y(), mc_rbdyn::rpyFromMat(pose.rotation()).z());
  };
  auto convertTo3d = [](const Eigen::Vector3d & trans) -> sva::PTransformd {
    return sva::PTransformd(sva::RotZ(trans.z()), Eigen::Vector3d(trans.x(), trans.y(), 0));
  };

  // The 2D variables (i.e., targetTrans, deltaTrans) represent the transformation relative to the initial pose,
  // while the 3D variables (i.e., initialFootMidpose, goalFootMidpose, footMidpose) represent the
  // transformation in the world frame.
  const sva::PTransformd & initialFootMidpose =
      projGround(sva::interpolate(targetFootPoses_.at(Foot::Left), targetFootPoses_.at(Foot::Right), 0.5));

  Foot foot = (waypointTransList.empty() ? targetTrans.y() : waypointTransList[0].y()) >= 0 ? Foot::Left : Foot::Right;
  sva::PTransformd footMidpose = initialFootMidpose;
  double startTime = ctl().t() + 1.0;

  for(size_t i = 0; i < waypointTransList.size() + 1; i++)
  {
    const Eigen::Vector3d & goalTrans = (i == waypointTransList.size() ? targetTrans : waypointTransList[i]);
    const sva::PTransformd & goalFootMidpose = convertTo3d(goalTrans) * initialFootMidpose;
    double thre = (i == waypointTransList.size() ? 1e-6 : 1e-2);

    while(convertTo2d(goalFootMidpose * footMidpose.inv()).norm() > thre)
    {
      Eigen::Vector3d deltaTrans = convertTo2d(goalFootMidpose * footMidpose.inv());
      footMidpose = convertTo3d(clampDeltaTrans(deltaTrans, foot)) * footMidpose;

      const auto & footstep = makeFootstep(foot, footMidpose, startTime);
      appendFootstep(footstep);

      foot = opposite(foot);
      startTime = footstep.transitEndTime;
    }
  }

  for(int i = 0; i < lastFootstepNum + 1; i++)
  {
    const auto & footstep = makeFootstep(foot, footMidpose, startTime);
    appendFootstep(footstep);

    foot = opposite(foot);
    startTime = footstep.transitEndTime;
  }

  return true;
}

bool FootManager::startVelMode()
{
  if(velModeData_.enabled_)
  {
    mc_rtc::log::warning("[FootManager] It is already in velocity mode, but startVelMode is called.");
    return false;
  }

  if(footstepQueue_.size() > 0)
  {
    mc_rtc::log::error("[FootManager] startVelMode is available only when the footstep queue is empty: {}",
                       footstepQueue_.size());
    return false;
  }

  velModeData_.reset(true);

  // Add footsteps to queue for walking in place
  Foot foot = Foot::Left;
  const sva::PTransformd & footMidpose =
      projGround(sva::interpolate(targetFootPoses_.at(Foot::Left), targetFootPoses_.at(Foot::Right), 0.5));
  double startTime = ctl().t() + 1.0;
  for(int i = 0; i < velModeData_.config_.footstepQueueSize; i++)
  {
    const auto & footstep = makeFootstep(foot, footMidpose, startTime);
    appendFootstep(footstep);

    foot = opposite(foot);
    startTime = footstep.transitEndTime;
  }

  return true;
}

bool FootManager::endVelMode()
{
  if(!velModeData_.enabled_)
  {
    mc_rtc::log::warning("[FootManager] It is not in velocity mode, but endVelMode is called.");
    return false;
  }

  velModeData_.reset(false);

  // Update last footstep pose to align both feet
  // Note that this process assumes that velModeData_.config_.footstepQueueSize is at least 3
  const auto & lastFootstep1 = *(footstepQueue_.rbegin() + 1);
  auto & lastFootstep2 = *(footstepQueue_.rbegin());
  sva::PTransformd footMidpose = config_.midToFootTranss.at(lastFootstep1.foot).inv() * lastFootstep1.pose;
  lastFootstep2.pose = config_.midToFootTranss.at(lastFootstep2.foot) * footMidpose;

  return true;
}
void FootManager::updateFootTraj()
{
  // デフォルトでホールドモードをオフ（足タスクの目標ジャンプ抑制）
  for(const auto & foot : Feet::Both)
  {
    ctl().footTasks_.at(foot)->hold(false);
  }

  // 古いフットステップ（既に終了したもの）をキューから削除
  while(!footstepQueue_.empty() && footstepQueue_.front().transitEndTime < ctl().t())
  {
    prevFootstep_ = std::make_shared<Footstep>(footstepQueue_.front());
    footstepQueue_.pop_front();
  }

  // スイングフェーズ中（現在時刻がスイング開始〜終了時間内）
  if(!footstepQueue_.empty() && footstepQueue_.front().swingStartTime <= ctl().t()
     && ctl().t() <= footstepQueue_.front().swingEndTime)
  {
    // シングルサポートフェーズ

    if(swingFootstep_)
    {
      // 既にスイングフットステップを持っている場合、整合性チェック
      if(swingFootstep_ != &(footstepQueue_.front()))
      {
        mc_rtc::log::error_and_throw("[FootManager] Swing footstep is not consistent.");
      }

      // 終端姿勢を同期（目標更新）
      swingFootstep_->pose = swingTraj_->endPose_;
    }
    else
    {
      // 新しくスイングフットステップをセット
      swingFootstep_ = &(footstepQueue_.front());

      // IK目標のジャンプを防ぐためホールドモードON
      ctl().footTasks_.at(swingFootstep_->foot)->hold(true);

      // スイング軌道生成
      {
        const sva::PTransformd & swingStartPose = ctl().robot().surfacePose(surfaceName(swingFootstep_->foot));
        sva::PTransformd swingEndPose = swingFootstep_->pose;

        // 終端ポーズ補正
        auto &p = swingEndPose.translation();

        // Z 方向（上正）の補正
        p.z() += landingZOffset_.at(swingFootstep_->foot);

        // 上書き着地姿勢指定時は、前足ステップを参照して終端姿勢を補正
        if(config_.overwriteLandingPose && prevFootstep_)
        {
          sva::PTransformd swingRelPose = swingFootstep_->pose * prevFootstep_->pose.inv();
          swingEndPose.translation() = (swingRelPose * targetFootPoses_.at(prevFootstep_->foot)).translation();
        }

        // --- apply landing offsets to end pose ---
        {
          auto & t = swingEndPose.translation();
          // Z
          t.z() += landingZOffset_.at(swingFootstep_->foot);
        }


        // スイング軌道タイプに応じて生成
        std::string swingTrajType =
            swingFootstep_->swingTrajConfig("type", static_cast<std::string>(config_.defaultSwingTrajType));
        if(swingTrajType == "CubicSplineSimple")
        {
          swingTraj_ = std::make_shared<SwingTrajCubicSplineSimple>(swingStartPose, swingEndPose,
              swingFootstep_->swingStartTime, swingFootstep_->swingEndTime, config_.footTaskGain, swingFootstep_->swingTrajConfig);
        }
        else if(swingTrajType == "IndHorizontalVertical")
        {
          swingFootstep_->swingTrajConfig.add("localVertexList",
              calcSurfaceVertexList(ctl().robot().surface(surfaceName(swingFootstep_->foot)), sva::PTransformd::Identity()));
          swingTraj_ = std::make_shared<SwingTrajIndHorizontalVertical>(swingStartPose, swingEndPose,
              swingFootstep_->swingStartTime, swingFootstep_->swingEndTime, config_.footTaskGain, swingFootstep_->swingTrajConfig);
        }
        else if(swingTrajType == "VariableTaskGain")
        {
          swingTraj_ = std::make_shared<SwingTrajVariableTaskGain>(swingStartPose, swingEndPose,
              swingFootstep_->swingStartTime, swingFootstep_->swingEndTime, config_.footTaskGain, swingFootstep_->swingTrajConfig);
        }
        else if(swingTrajType == "LandingSearch")
        {
          swingTraj_ = std::make_shared<SwingTrajLandingSearch>(swingStartPose, swingEndPose,
              swingFootstep_->swingStartTime, swingFootstep_->swingEndTime, config_.footTaskGain, swingFootstep_->swingTrajConfig);
        }
        else
        {
          mc_rtc::log::error_and_throw("[FootManager] Invalid swingTrajType: {}.", swingTrajType);
        }
      }

      // baseYaw補間関数を生成
      {
        double swingStartBaseYaw =
            mc_rbdyn::rpyFromMat(
                TrajColl::interpolate<Eigen::Matrix3d>(
                    targetFootPoses_.at(Foot::Left).rotation().transpose(),
                    targetFootPoses_.at(Foot::Right).rotation().transpose(), 0.5).transpose()).z();
        baseYawFunc_->appendPoint(std::make_pair(swingFootstep_->swingStartTime,
            Eigen::AngleAxisd(swingStartBaseYaw, Eigen::Vector3d::UnitZ()).toRotationMatrix()));

        double swingEndBaseYaw =
            mc_rbdyn::rpyFromMat(
                TrajColl::interpolate<Eigen::Matrix3d>(
                    swingFootstep_->pose.rotation().transpose(),
                    targetFootPoses_.at(opposite(swingFootstep_->foot)).rotation().transpose(), 0.5).transpose()).z();
        baseYawFunc_->appendPoint(std::make_pair(swingFootstep_->swingEndTime,
            Eigen::AngleAxisd(swingEndBaseYaw, Eigen::Vector3d::UnitZ()).toRotationMatrix()));

        baseYawFunc_->calcCoeff();
      }

      // アームスイング補間器設定（必要に応じ）
      if(config_.enableArmSwing)
      {
        // 条件によりアームスイング生成
        int totalSize = 0;
        for(const auto & jointAngleKV : config_.jointAnglesForArmSwing.at("Nominal"))
        {
          totalSize += static_cast<int>(jointAngleKV.second.size());
        }
        sva::PTransformd startToEndTrans = swingTraj_->endPose_ * swingTraj_->startPose_.inv();
        double forwardDist = startToEndTrans.translation().x();
        double forwardAngle = std::abs(std::atan2(startToEndTrans.translation().y(), startToEndTrans.translation().x()));
        constexpr double forwardDistThre = 0.1; // [m]
        constexpr double forwardAngleThre = mc_rtc::constants::toRad(30.0); // [rad]
        if(totalSize > 0 && forwardDist > forwardDistThre && forwardAngle < forwardAngleThre)
        {
          // CubicSpline補間器へアーム関節角補間点追加
          auto jointAnglesMapToVec = [totalSize](const std::map<std::string, std::vector<double>> & jointAnglesMap) -> Eigen::VectorXd {
            Eigen::VectorXd jointAnglesVec(totalSize);
            int vecIdx = 0;
            for(const auto & jointAngleKV : jointAnglesMap)
            {
              jointAnglesVec.segment(vecIdx, jointAngleKV.second.size()) =
                  Eigen::Map<const Eigen::VectorXd>(jointAngleKV.second.data(), jointAngleKV.second.size());
              vecIdx += static_cast<int>(jointAngleKV.second.size());
            }
            return jointAnglesVec;
          };
          TrajColl::BoundaryConstraint<Eigen::VectorXd> zeroVelBC(TrajColl::BoundaryConstraintType::Velocity, Eigen::VectorXd::Zero(totalSize));
          std::map<std::string, std::vector<double>> currentJointAnglesMap;
          auto postureTask = ctl().getPostureTask(ctl().robot().name());
          for(const auto & jointAngleKV : config_.jointAnglesForArmSwing.at("Nominal"))
          {
            currentJointAnglesMap[jointAngleKV.first] =
                postureTask->posture()[ctl().robot().jointIndexByName(jointAngleKV.first)];
          }
          Eigen::VectorXd currentJointAnglesVec = jointAnglesMapToVec(currentJointAnglesMap);
          Eigen::VectorXd swingJointAnglesVec = jointAnglesMapToVec(config_.jointAnglesForArmSwing.at(std::to_string(swingFootstep_->foot)));
          Eigen::VectorXd nominalJointAnglesVec = jointAnglesMapToVec(config_.jointAnglesForArmSwing.at("Nominal"));
          armSwingFunc_ = std::make_shared<TrajColl::CubicSpline<Eigen::VectorXd>>(totalSize, zeroVelBC, zeroVelBC);
          armSwingFunc_->appendPoint(std::make_pair(swingFootstep_->swingStartTime, currentJointAnglesVec));
          armSwingFunc_->appendPoint(std::make_pair(0.5 * (swingFootstep_->swingStartTime + swingFootstep_->swingEndTime), swingJointAnglesVec));
          armSwingFunc_->appendPoint(std::make_pair(swingFootstep_->swingEndTime, 0.5 * (nominalJointAnglesVec + swingJointAnglesVec)));
          armSwingFunc_->appendPoint(std::make_pair(swingFootstep_->transitEndTime, nominalJointAnglesVec));
          armSwingFunc_->calcCoeff();
        }
      }

      // サポートフェーズ状態更新
      supportPhase_ = (swingFootstep_->foot == Foot::Left) ? SupportPhase::RightSupport : SupportPhase::LeftSupport;
    }

    // タッチダウン判定
    if(!touchDown_ && detectTouchDown())
    {
      touchDown_ = true;

      if(config_.stopSwingTrajForTouchDownFoot)
      {
        swingTraj_->touchDown(ctl().t());
      }
      // --- added: clear offset of the touchdown foot ---
      this->clearLandingZOffset(swingFootstep_->foot);
    }

    // スイング軌道から現在のターゲット姿勢・速度・加速度・ゲインを取得して設定
    {
      swingTraj_->update(ctl().t());
      targetFootPoses_.at(swingFootstep_->foot) = swingTraj_->pose(ctl().t());
      targetFootVels_.at(swingFootstep_->foot) = swingTraj_->vel(ctl().t());
      targetFootAccels_.at(swingFootstep_->foot) = swingTraj_->accel(ctl().t());
      footTaskGains_.at(swingFootstep_->foot) = swingTraj_->taskGain(ctl().t());
    }
  }
  else
  {
    // ダブルサポートフェーズ処理
    if(swingFootstep_)
    {
      if(!(config_.keepPoseForTouchDownFoot && touchDown_))
      {
        targetFootPoses_.at(swingFootstep_->foot) = swingTraj_->endPose_;
      }
      targetFootVels_.at(swingFootstep_->foot) = sva::MotionVecd::Zero();
      targetFootAccels_.at(swingFootstep_->foot) = sva::MotionVecd::Zero();
      footTaskGains_.at(swingFootstep_->foot) = config_.footTaskGain;

      // ★ここから追加：DS中に世界基準へzを戻す補間を両足に設定
      {
        const double t0 = ctl().t();
        const double t1 = swingFootstep_->transitEndTime; // DS 終了時刻

        // 現在ターゲットの左右 z と平均
        double zL0 = targetFootPoses_.at(Foot::Left).translation().z();
        double zR0 = targetFootPoses_.at(Foot::Right).translation().z();
        double curMidZ = 0.5 * (zL0 + zR0);

        // 世界基準との差分を左右に0.5ずつ配分（1歩あたりの戻し量に上限）
        double dz = worldGroundZ0_ - curMidZ;
        const double dz_limit = 0.015; // 例：±15mm/歩 に制限
        double dz_each = std::clamp(0.5 * dz, -dz_limit, dz_limit);

        // 左右それぞれ終点 z
        double zL1 = zL0 + dz_each;
        double zR1 = zR0 + dz_each;

        // スイング足: 既存の補間器を作る部分を置き換え（終点zだけ z?1 に）
        {
          auto pose0 = swingTraj_->endPose_;
          auto pose1 = targetFootPoses_.at(swingFootstep_->foot);
          // 現在のターゲットをベースに z 終点を上書き
          auto p1 = pose1.translation();
          p1.z() = (swingFootstep_->foot == Foot::Left) ? zL1 : zR1;
          pose1.translation() = p1;

          auto trajStartFootPoseFunc =
            std::make_shared<TrajColl::CubicInterpolator<sva::PTransformd, sva::MotionVecd>>();
          trajStartFootPoseFunc->appendPoint(std::make_pair(t0, pose0));
          trajStartFootPoseFunc->appendPoint(std::make_pair(t1, pose1));
          trajStartFootPoseFunc->calcCoeff();
          trajStartFootPoseFuncs_.at(swingFootstep_->foot) = trajStartFootPoseFunc;

          // DS中のターゲットも一貫させる
          targetFootPoses_.at(swingFootstep_->foot) = pose1;
        }

        // 支持足側にも z 復帰用の補間器を設定
        {
          Foot support = opposite(swingFootstep_->foot);
          auto pose0 = targetFootPoses_.at(support);
          auto pose1 = pose0;
          auto p = pose1.translation();
          p.z() = (support == Foot::Left) ? zL1 : zR1;
          pose1.translation() = p;

          auto trajStartFootPoseFunc =
            std::make_shared<TrajColl::CubicInterpolator<sva::PTransformd, sva::MotionVecd>>();
          trajStartFootPoseFunc->appendPoint(std::make_pair(t0, pose0));
          trajStartFootPoseFunc->appendPoint(std::make_pair(t1, pose1));
          trajStartFootPoseFunc->calcCoeff();
          trajStartFootPoseFuncs_.at(support) = trajStartFootPoseFunc;

          targetFootPoses_.at(support) = pose1;
        }
      }
      // ★ここまで追加

      supportPhase_ = SupportPhase::DoubleSupport;
      swingTraj_.reset();
      baseYawFunc_->clearPoints();
      touchDown_ = false;
      swingFootstep_ = nullptr;
    }
  }


  // 足タスクへ現在のターゲットをセット
  for(const auto & foot : Feet::Both)
  {
    ctl().footTasks_.at(foot)->targetPose(targetFootPoses_.at(foot));
    ctl().footTasks_.at(foot)->targetVel(targetFootVels_.at(foot));
    ctl().footTasks_.at(foot)->targetAccel(targetFootAccels_.at(foot));
    ctl().footTasks_.at(foot)->setGains(footTaskGains_.at(foot).stiffness, footTaskGains_.at(foot).damping);
  }

  // インピーダンスゲイン種別の更新
  {
    std::unordered_map<Foot, std::string> newImpGainTypes;
    const auto & contactFeet = getCurrentContactFeet();
    if(contactFeet.size() == 1)
    {
      newImpGainTypes.emplace(*(contactFeet.cbegin()), "SingleSupport");
      newImpGainTypes.emplace(opposite(*(contactFeet.cbegin())), "Swing");
    }
    else
    {
      for(const auto & foot : Feet::Both)
      {
        newImpGainTypes.emplace(foot, "DoubleSupport");
      }
    }
    for(const auto & foot : Feet::Both)
    {
      if(requireImpGainUpdate_)
      {
        break;
      }
      if(impGainTypes_.at(foot) != newImpGainTypes.at(foot))
      {
        requireImpGainUpdate_ = true;
      }
    }
    impGainTypes_ = newImpGainTypes;
  }

  // インピーダンスゲイン適用
  if(requireImpGainUpdate_)
  {
    requireImpGainUpdate_ = false;
    for(const auto & foot : Feet::Both)
    {
      ctl().footTasks_.at(foot)->gains() = config_.impGains.at(impGainTypes_.at(foot));
    }
  }

  // ベース姿勢タスク更新（フェーズに応じて）
  if(supportPhase_ == SupportPhase::DoubleSupport)
  {
    const sva::PTransformd & footMidpose = sva::interpolate(targetFootPoses_.at(Foot::Left), targetFootPoses_.at(Foot::Right), 0.5);
    ctl().baseOriTask_->orientation(sva::RotZ(mc_rbdyn::rpyFromMat(footMidpose.rotation()).z()));
    ctl().baseOriTask_->refVel(Eigen::Vector3d::Zero());
    ctl().baseOriTask_->refAccel(Eigen::Vector3d::Zero());
  }
  else
  {
    ctl().baseOriTask_->orientation((*baseYawFunc_)(ctl().t()).transpose());
    ctl().baseOriTask_->refVel(baseYawFunc_->derivative(ctl().t(), 1));
    ctl().baseOriTask_->refAccel(baseYawFunc_->derivative(ctl().t(), 2));
  }

  // アームスイング更新
  if(armSwingFunc_)
  {
    if(armSwingFunc_->domainUpperLimit() < ctl().t())
    {
      armSwingFunc_.reset();
    }
    else
    {
      auto jointAnglesVecToMap = [this](const Eigen::VectorXd & jointAnglesVec) -> std::map<std::string, std::vector<double>> {
        std::map<std::string, std::vector<double>> jointAnglesMap;
        int vecIdx = 0;
        for(const auto & jointAngleKV : config_.jointAnglesForArmSwing.at("Nominal"))
        {
          jointAnglesMap[jointAngleKV.first] = std::vector<double>(
              jointAnglesVec.data() + vecIdx, jointAnglesVec.data() + vecIdx + jointAngleKV.second.size());
          vecIdx += static_cast<int>(jointAngleKV.second.size());
        }
        return jointAnglesMap;
      };
      auto postureTask = ctl().getPostureTask(ctl().robot().name());
      postureTask->target(jointAnglesVecToMap((*armSwingFunc_)(ctl().t())));
    }
  }

  // GUI用にフットステップ表示を更新
  {
    std::vector<std::vector<Eigen::Vector3d>> footstepPolygonList;
    for(const auto & footstep : footstepQueue_)
    {
      const auto & surface = ctl().robot().surface(surfaceName(footstep.foot));
      footstepPolygonList.push_back(calcSurfaceVertexList(surface, footstep.pose));
    }
    ctl().gui()->removeCategory({ctl().name(), config_.name, "FootstepMarker"});
    ctl().gui()->addElement({ctl().name(), config_.name, "FootstepMarker"},
                            mc_rtc::gui::Polygon("Footstep", {mc_rtc::gui::Color::Blue, 0.02},
                                                 [footstepPolygonList]() { return footstepPolygonList; }));
  }
}

void FootManager::updateZmpTraj()
{
  zmpFunc_->clearPoints();
  groundPosZFunc_->clearPoints();
  contactFootPosesList_.clear();

  auto groundZ = [this](const std::unordered_map<Foot, sva::PTransformd> & _footPoses) {
    if(keepWorldGroundBaseline_) { return worldGroundZ0_; }
    // 従来の足裏平均を使いたい場合
    return 0.5 * (_footPoses.at(Foot::Left).translation().z()
                + _footPoses.at(Foot::Right).translation().z());
  };

  // Update trajStartFootPoses_
  for(auto & trajStartFootPoseFuncKV : trajStartFootPoseFuncs_)
  {
    auto & trajStartFootPoseFunc = trajStartFootPoseFuncKV.second;
    if(!trajStartFootPoseFunc)
    {
      continue;
    }
    trajStartFootPoses_.at(trajStartFootPoseFuncKV.first) =
        (*trajStartFootPoseFunc)(std::min(ctl().t(), trajStartFootPoseFunc->endTime()));
    if(trajStartFootPoseFunc->endTime() <= ctl().t())
    {
      trajStartFootPoseFunc.reset();
    }
  }
  std::unordered_map<Foot, sva::PTransformd> footPoses = trajStartFootPoses_;

  auto calcFootMidposZ = [](const std::unordered_map<Foot, sva::PTransformd> & _footPoses) {
    return 0.5 * (_footPoses.at(Foot::Left).translation().z() + _footPoses.at(Foot::Right).translation().z());
  };

  if(footstepQueue_.empty() || ctl().t() < footstepQueue_.front().transitStartTime)
  {
    // Set initial point
    zmpFunc_->appendPoint(std::make_pair(ctl().t(), calcZmpWithOffset(footPoses)));
    groundPosZFunc_->appendPoint(std::make_pair(ctl().t(), calcFootMidposZ(footPoses)));
    contactFootPosesList_.emplace(ctl().t(), footPoses);
  }

  for(const auto & footstep : footstepQueue_)
  {
    Foot supportFoot = opposite(footstep.foot);
    Eigen::Vector3d supportFootZmp = calcZmpWithOffset(supportFoot, footPoses.at(supportFoot));

    if(ctl().t() <= footstep.swingEndTime)
    {
      // zmpFunc_->appendPoint(std::make_pair(footstep.transitStartTime, calcZmpWithOffset(footPoses)));
      // groundPosZFunc_->appendPoint(std::make_pair(footstep.transitStartTime, calcFootMidposZ(footPoses)));
      // contactFootPosesList_.emplace(footstep.transitStartTime, footPoses);

      // zmpFunc_->appendPoint(std::make_pair(footstep.swingStartTime, supportFootZmp));
      // groundPosZFunc_->appendPoint(std::make_pair(footstep.swingStartTime, calcFootMidposZ(footPoses)));
      // contactFootPosesList_.emplace(footstep.swingStartTime, std::unordered_map<Foot, sva::PTransformd>{
      //                                                            {supportFoot, footPoses.at(supportFoot)}});

      zmpFunc_->appendPoint(std::make_pair(footstep.transitStartTime, calcZmpWithOffset(footPoses)));
      groundPosZFunc_->appendPoint(std::make_pair(footstep.transitStartTime, groundZ(footPoses)));
      contactFootPosesList_.emplace(footstep.transitStartTime, footPoses);

      zmpFunc_->appendPoint(std::make_pair(footstep.swingStartTime, supportFootZmp));
      groundPosZFunc_->appendPoint(std::make_pair(footstep.swingStartTime, groundZ(footPoses)));
      contactFootPosesList_.emplace(footstep.swingStartTime,
        std::unordered_map<Foot, sva::PTransformd>{{supportFoot, footPoses.at(supportFoot)}});
      

      // Update footPoses
      footPoses.at(footstep.foot) = (footstep.swingStartTime <= ctl().t() ? swingTraj_->endPose_ : footstep.pose);
    }

    // zmpFunc_->appendPoint(std::make_pair(footstep.swingEndTime, supportFootZmp));
    // groundPosZFunc_->appendPoint(std::make_pair(footstep.swingEndTime, calcFootMidposZ(footPoses)));
    // contactFootPosesList_.emplace(footstep.swingEndTime, footPoses);

    // groundPosZFunc_->appendPoint(std::make_pair(footstep.transitEndTime, calcFootMidposZ(footPoses)));
    // zmpFunc_->appendPoint(std::make_pair(footstep.transitEndTime, calcZmpWithOffset(footPoses)));
    // contactFootPosesList_.emplace(footstep.transitEndTime, footPoses);

    zmpFunc_->appendPoint(std::make_pair(footstep.swingEndTime, supportFootZmp));
    groundPosZFunc_->appendPoint(std::make_pair(footstep.swingEndTime, groundZ(footPoses)));
    contactFootPosesList_.emplace(footstep.swingEndTime, footPoses);

    groundPosZFunc_->appendPoint(std::make_pair(footstep.transitEndTime, groundZ(footPoses)));
    zmpFunc_->appendPoint(std::make_pair(footstep.transitEndTime, calcZmpWithOffset(footPoses)));
    contactFootPosesList_.emplace(footstep.transitEndTime, footPoses);

    if(ctl().t() + config_.zmpHorizon <= footstep.transitEndTime)
    {
      break;
    }
  }

  // if(footstepQueue_.empty() || footstepQueue_.back().transitEndTime < ctl().t() + config_.zmpHorizon)
  // {
  //   // Set terminal point
  //   zmpFunc_->appendPoint(std::make_pair(ctl().t() + config_.zmpHorizon, calcZmpWithOffset(footPoses)));
  //   groundPosZFunc_->appendPoint(std::make_pair(ctl().t() + config_.zmpHorizon, calcFootMidposZ(footPoses)));
  // }

  if(footstepQueue_.empty() || footstepQueue_.back().transitEndTime < ctl().t() + config_.zmpHorizon)
  {
    zmpFunc_->appendPoint(std::make_pair(ctl().t() + config_.zmpHorizon, calcZmpWithOffset(footPoses)));
    groundPosZFunc_->appendPoint(std::make_pair(ctl().t() + config_.zmpHorizon, groundZ(footPoses)));
  }

  zmpFunc_->calcCoeff();
  groundPosZFunc_->calcCoeff();
}

void FootManager::updateVelMode()
{
  auto convertTo2d = [](const sva::PTransformd & pose) -> Eigen::Vector3d {
    return Eigen::Vector3d(pose.translation().x(), pose.translation().y(), mc_rbdyn::rpyFromMat(pose.rotation()).z());
  };
  auto convertTo3d = [](const Eigen::Vector3d & trans) -> sva::PTransformd {
    return sva::PTransformd(sva::RotZ(trans.z()), Eigen::Vector3d(trans.x(), trans.y(), 0));
  };

  // Keep the next footstep and delete the second and subsequent footsteps
  footstepQueue_.erase(footstepQueue_.begin() + 1, footstepQueue_.end());
  const auto & nextFootstep = footstepQueue_.front();
  sva::PTransformd footMidpose = projGround(config_.midToFootTranss.at(nextFootstep.foot).inv() * nextFootstep.pose);
  Eigen::Vector3d deltaTrans = config_.footstepDuration * velModeData_.targetVel_;

  // Update footstep online during swing
  if(velModeData_.config_.enableOnlineFootstepUpdate && swingTraj_ && swingTraj_->type() == "VariableTaskGain")
  {
    constexpr double updateEndTimeRatio = 0.9;
    double approachTime = std::dynamic_pointer_cast<SwingTrajVariableTaskGain>(swingTraj_)->approachTime_;
    double updateEndTime = updateEndTimeRatio * (approachTime - swingTraj_->startTime_) + swingTraj_->startTime_;
    if(swingTraj_->startTime_ <= ctl().t() && ctl().t() <= updateEndTime)
    {
      sva::PTransformd currentFootMidpose = projGround(config_.midToFootTranss.at(opposite(nextFootstep.foot)).inv()
                                                       * targetFootPoses_.at(opposite(nextFootstep.foot)));
      footMidpose = convertTo3d(clampDeltaTrans(deltaTrans, nextFootstep.foot)) * currentFootMidpose;
      const sva::PTransformd & footstepPoseOrig = nextFootstep.pose;
      sva::PTransformd footstepPoseNew = config_.midToFootTranss.at(nextFootstep.foot) * footMidpose;
      Eigen::Vector3d footstepTrans = convertTo2d(footstepPoseNew * footstepPoseOrig.inv());
      double remainingDurationRatio =
          std::clamp((updateEndTime - ctl().t()) / (updateEndTime - swingTraj_->startTime_), 0.0, 1.0);
      Eigen::Vector3d footstepTransMax =
          std::pow(remainingDurationRatio, 2) * Eigen::Vector3d(0.1, 0.1, mc_rtc::constants::toRad(15));
      Eigen::Vector3d footstepTransMin = -1 * footstepTransMax;
      sva::PTransformd footstepPoseNewClamped =
          convertTo3d(mc_filter::utils::clamp(footstepTrans, footstepTransMin, footstepTransMax)) * footstepPoseOrig;
      swingTraj_->endPose_ = footstepPoseNewClamped;
    }
  }

  // Append new footsteps
  Foot foot = opposite(nextFootstep.foot);
  double startTime = nextFootstep.transitEndTime;
  for(int i = 0; i < velModeData_.config_.footstepQueueSize - 1; i++)
  {
    footMidpose = convertTo3d(clampDeltaTrans(deltaTrans, foot)) * footMidpose;

    const auto & footstep = makeFootstep(foot, footMidpose, startTime);
    footstepQueue_.push_back(footstep);

    foot = opposite(foot);
    startTime = footstep.transitEndTime;
  }
}

double FootManager::touchDownRemainingDuration() const
{
  if(supportPhase_ == SupportPhase::DoubleSupport)
  {
    return 0;
  }
  else
  {
    return swingFootstep_->swingEndTime - ctl().t();
  }
}

bool FootManager::detectTouchDown() const
{
  // False for double support phase
  if(supportPhase_ == SupportPhase::DoubleSupport)
  {
    return false;
  }

  // False if the remaining duration does not meet the threshold
  if(touchDownRemainingDuration() > config_.touchDownRemainingDuration)
  {
    return false;
  }

  // False if the position error does not meet the threshold
  if((swingTraj_->endPose_.translation() - swingTraj_->pose(ctl().t()).translation()).norm()
     > config_.touchDownPosError)
  {
    return false;
  }

  // False if the normal force does not meet the threshold
  Foot swingFoot = (supportPhase_ == SupportPhase::LeftSupport ? Foot::Right : Foot::Left);
  double fz = ctl().robot().surfaceWrench(surfaceName(swingFoot)).force().z();
  if(fz < config_.touchDownForceZ)
  {
    return false;
  }

  return true;
}
