#ifdef __CINT__

#pragma link C++ class art::EventAuxiliary+;
#pragma link C++ class art::EventID+;
#pragma link C++ class art::SubRunID+;
#pragma link C++ class art::RunID+;
#pragma link C++ class art::Timestamp+;
#pragma link C++ class art::SubRunAuxiliary+;
#pragma link C++ class art::Hash<2>+;
#pragma link C++ class art::RunAuxiliary+;
#pragma link C++ class art::ResultsAuxiliary+;
#pragma link C++ class art::History+;
#pragma link C++ class art::Wrapper<art::TriggerResults>+;
#pragma link C++ class art::EDProduct+;
#pragma link C++ class art::TriggerResults+;
#pragma link C++ class art::HLTGlobalStatus+;
#pragma link C++ class fhicl::ParameterSetID+;
#pragma link C++ class art::Wrapper<vector<artdaq::Fragment> >+;
#pragma link C++ class artdaq::Fragment+;
#pragma link C++ class artdaq::QuickVec<ULong64_t>+;
#pragma link C++ class art::Wrapper<vector<artdaq::PackageBuildInfo> >+;
#pragma link C++ class artdaq::PackageBuildInfo+;
#pragma link C++ class art::ProductProvenance+;
#pragma link C++ class art::ProductID+;
#pragma link C++ class art::Hash<5>+;
#pragma link C++ class art::Transient<art::ProductProvenance::Transients>+;
#pragma link C++ class art::HLTPathStatus+;
#pragma link C++ class art::FileFormatVersion+;
#pragma link C++ class art::FileIndex::Element+;
#pragma link C++ class art::ProcessHistory+;
#pragma link C++ class art::ProcessConfiguration+;
#pragma link C++ class art::Transient<art::ProcessHistory::Transients>+;
#pragma link C++ class art::ProductRegistry+;
#pragma link C++ class art::BranchKey+;
#pragma link C++ class art::BranchDescription+;
#pragma link C++ class art::Hash<3>+;
#pragma link C++ class art::Transient<art::BranchDescription::Transients>+;
#pragma link C++ class art::Parentage+;
#pragma link C++ class art::BranchChildren+;
#pragma link C++ class art::Wrapper<art::Assns<CRT::Trigger,anab::CosmicTag,void> >+;
#pragma link C++ class art::Assns<CRT::Trigger,anab::CosmicTag,void>+;
#pragma link C++ class art::detail::AssnsBase+;
#pragma link C++ class art::Wrapper<art::Assns<anab::CosmicTag,anab::T0,void> >+;
#pragma link C++ class art::Assns<anab::CosmicTag,anab::T0,void>+;
#pragma link C++ class art::Wrapper<art::Assns<raw::RawDigit,raw::RDTimeStamp,void> >+;
#pragma link C++ class art::Assns<raw::RawDigit,raw::RDTimeStamp,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Cluster,recob::Hit,void> >+;
#pragma link C++ class art::Assns<recob::Cluster,recob::Hit,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Hit,recob::SpacePoint,void> >+;
#pragma link C++ class art::Assns<recob::Hit,recob::SpacePoint,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::OpFlash,recob::OpHit,void> >+;
#pragma link C++ class art::Assns<recob::OpFlash,recob::OpHit,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::PFParticle,anab::T0,void> >+;
#pragma link C++ class art::Assns<recob::PFParticle,anab::T0,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::PFParticle,larpandoraobj::PFParticleMetadata,void> >+;
#pragma link C++ class art::Assns<recob::PFParticle,larpandoraobj::PFParticleMetadata,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::PFParticle,recob::Cluster,void> >+;
#pragma link C++ class art::Assns<recob::PFParticle,recob::Cluster,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::PFParticle,recob::PCAxis,void> >+;
#pragma link C++ class art::Assns<recob::PFParticle,recob::PCAxis,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::PFParticle,recob::Shower,void> >+;
#pragma link C++ class art::Assns<recob::PFParticle,recob::Shower,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::PFParticle,recob::Slice,void> >+;
#pragma link C++ class art::Assns<recob::PFParticle,recob::Slice,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::PFParticle,recob::SpacePoint,void> >+;
#pragma link C++ class art::Assns<recob::PFParticle,recob::SpacePoint,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::PFParticle,recob::Track,void> >+;
#pragma link C++ class art::Assns<recob::PFParticle,recob::Track,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::PFParticle,recob::Vertex,void> >+;
#pragma link C++ class art::Assns<recob::PFParticle,recob::Vertex,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Shower,anab::Calorimetry,void> >+;
#pragma link C++ class art::Assns<recob::Shower,anab::Calorimetry,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Shower,recob::Hit,void> >+;
#pragma link C++ class art::Assns<recob::Shower,recob::Hit,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Shower,recob::PCAxis,void> >+;
#pragma link C++ class art::Assns<recob::Shower,recob::PCAxis,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Slice,recob::Hit,void> >+;
#pragma link C++ class art::Assns<recob::Slice,recob::Hit,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::SpacePoint,recob::Hit,void> >+;
#pragma link C++ class art::Assns<recob::SpacePoint,recob::Hit,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Track,anab::Calorimetry,void> >+;
#pragma link C++ class art::Assns<recob::Track,anab::Calorimetry,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Track,anab::CosmicTag,void> >+;
#pragma link C++ class art::Assns<recob::Track,anab::CosmicTag,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Track,anab::ParticleID,void> >+;
#pragma link C++ class art::Assns<recob::Track,anab::ParticleID,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Track,anab::T0,void> >+;
#pragma link C++ class art::Assns<recob::Track,anab::T0,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Track,recob::Hit,recob::TrackHitMeta> >+;
#pragma link C++ class art::Assns<recob::Track,recob::Hit,recob::TrackHitMeta>+;
#pragma link C++ class art::Assns<recob::Track,recob::Hit,void>+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Track,recob::Hit,void> >+;
#pragma link C++ class art::Wrapper<art::Assns<recob::Wire,recob::Hit,void> >+;
#pragma link C++ class art::Assns<recob::Wire,recob::Hit,void>+;
#pragma link C++ class art::Wrapper<vector<CRT::Trigger> >+;
#pragma link C++ class CRT::Trigger+;
#pragma link C++ class art::Wrapper<vector<anab::Calorimetry> >+;
#pragma link C++ class anab::Calorimetry+;
#pragma link C++ class geo::PlaneID+;
#pragma link C++ class geo::TPCID+;
#pragma link C++ class geo::CryostatID+;
#pragma link C++ class art::Wrapper<vector<anab::CosmicTag> >+;
#pragma link C++ class anab::CosmicTag+;
#pragma link C++ class art::Wrapper<vector<anab::FeatureVector<4> > >+;
#pragma link C++ class anab::FeatureVector<4>+;
#pragma link C++ class art::Wrapper<vector<anab::MVADescription<4> > >+;
#pragma link C++ class anab::MVADescription<4>+;
#pragma link C++ class art::Wrapper<vector<anab::ParticleID> >+;
#pragma link C++ class anab::ParticleID+;
#pragma link C++ class art::Wrapper<vector<anab::T0> >+;
#pragma link C++ class anab::T0+;
#pragma link C++ class art::Wrapper<vector<beam::ProtoDUNEBeamEvent> >+;
#pragma link C++ class beam::ProtoDUNEBeamEvent+;
#pragma link C++ class beam::CKov+;
#pragma link C++ class art::Wrapper<vector<dune::ProtoDUNETimeStamp> >+;
#pragma link C++ class dune::ProtoDUNETimeStamp+;
#pragma link C++ class art::Wrapper<vector<larpandoraobj::PFParticleMetadata> >+;
#pragma link C++ class larpandoraobj::PFParticleMetadata+;
#pragma link C++ class art::Wrapper<vector<raw::OpDetWaveform> >+;
#pragma link C++ class raw::OpDetWaveform+;
#pragma link C++ class art::Wrapper<vector<raw::RDStatus> >+;
#pragma link C++ class raw::RDStatus+;
#pragma link C++ class art::Wrapper<vector<raw::RDTimeStamp> >+;
#pragma link C++ class raw::RDTimeStamp+;
#pragma link C++ class art::Wrapper<vector<raw::ctb::pdspctb> >+;
#pragma link C++ class raw::ctb::pdspctb+;
#pragma link C++ class art::Wrapper<vector<recob::Cluster> >+;
#pragma link C++ class recob::Cluster+;
#pragma link C++ class art::Wrapper<vector<recob::Hit> >+;
#pragma link C++ class recob::Hit+;
#pragma link C++ class geo::WireID+;
#pragma link C++ class art::Wrapper<vector<recob::OpFlash> >+;
#pragma link C++ class recob::OpFlash+;
#pragma link C++ class art::Wrapper<vector<recob::OpHit> >+;
#pragma link C++ class recob::OpHit+;
#pragma link C++ class art::Wrapper<vector<recob::PCAxis> >+;
#pragma link C++ class recob::PCAxis+;
#pragma link C++ class art::Wrapper<vector<recob::PFParticle> >+;
#pragma link C++ class recob::PFParticle+;
#pragma link C++ class art::Wrapper<vector<recob::PointCharge> >+;
#pragma link C++ class recob::PointCharge+;
#pragma link C++ class art::Wrapper<vector<recob::Shower> >+;
#pragma link C++ class recob::Shower+;
#pragma link C++ class art::Wrapper<vector<recob::Slice> >+;
#pragma link C++ class recob::Slice+;
#pragma link C++ class art::Wrapper<vector<recob::SpacePoint> >+;
#pragma link C++ class recob::SpacePoint+;
#pragma link C++ class art::Wrapper<vector<recob::Track> >+;
#pragma link C++ class recob::Track+;
#pragma link C++ class recob::TrackTrajectory+;
#pragma link C++ class recob::Trajectory+;
#pragma link C++ class art::Wrapper<vector<recob::Vertex> >+;
#pragma link C++ class recob::Vertex+;
#pragma link C++ class art::Wrapper<vector<recob::Wire> >+;
#pragma link C++ class recob::Wire+;
#pragma link C++ class lar::sparse_vector<float>+;
#pragma link C++ class pair<art::RefCore,unsigned long>+;
#pragma link C++ class art::RefCore+;
#pragma link C++ class art::RefCore::RefCoreTransients+;
#pragma link C++ class recob::TrackHitMeta+;
#pragma link C++ class CRT::Hit+;
#pragma link C++ class beam::FBM+;
#pragma link C++ class recob::TrajectoryPointFlags+;
#pragma link C++ class util::flags::FlagSet<32,unsigned int>+;
#pragma link C++ class util::flags::BitMask<unsigned int>+;
#pragma link C++ class util::flags::Bits_t<unsigned int>+;
#pragma link C++ class raw::ctb::Trigger+;
#pragma link C++ class raw::ctb::ChStatus+;
#pragma link C++ class raw::ctb::Misc+;
#pragma link C++ class raw::ctb::WordIndex+;
//#pragma link C++ class lar::sparse_vector<float,::datarange_t>::datarange_t+;
#pragma link C++ class lar::range_t<unsigned long>+;
#pragma link C++ class art+;
//#pragma link C++ class art::EventAuxiliary::ExperimentType+;
#pragma link C++ class fhicl+;
#pragma link C++ class artdaq+;
//#pragma link C++ class art::ProductProvenance::Transients+;
#pragma link C++ class art::FileIndex+;
//#pragma link C++ class art::ProcessHistory::Transients+;
//#pragma link C++ class art::BranchType+;
//#pragma link C++ class art::BranchDescription::Transients+;
#pragma link C++ class CRT+;
#pragma link C++ class anab+;
#pragma link C++ class art::detail+;
#pragma link C++ class raw+;
//#pragma link C++ class raw::RawDigit+;
#pragma link C++ class recob+;
#pragma link C++ class larpandoraobj+;
#pragma link C++ class geo+;
//#pragma link C++ class anab::cosmic_tag_id+;
#pragma link C++ class beam+;
#pragma link C++ class dune+;
//#pragma link C++ class dune::ProtoDUNETimingCommand+;
#pragma link C++ class raw::ctb+;
//#pragma link C++ class raw::ctb::Feedback+;
//#pragma link C++ class geo::_plane_proj+;
//#pragma link C++ class geo::_plane_sigtype+;
//#pragma link C++ class recob::Vertex::Status+;
#pragma link C++ class lar+;
#pragma link C++ class util+;
#pragma link C++ class util::flags+;
#endif
