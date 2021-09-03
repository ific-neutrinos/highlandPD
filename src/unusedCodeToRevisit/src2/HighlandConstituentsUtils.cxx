#include "HighlandAnalysisUtils.hxx"
#include <stdio.h>
#include <math.h>


//**************************************************
int anaUtils::GetAllTracksUsingDet(const AnaBunchB& bunch, SubDetId::SubDetEnum det, AnaTrackB* selTracks[]) {
//**************************************************

  int count = 0;
  for (unsigned int it = 0; it < bunch.Particles.size(); it++) {
    AnaTrackB* track = static_cast<AnaTrackB*>(bunch.Particles[it]);
    if (anaUtils::TrackUsesDet(*track, det)) {
      selTracks[count] = track;
      count++;
    }
  }

  return count;
}


//**************************************************
int anaUtils::GetAllTracksUsingOnlyDet(const AnaBunchB& bunch, SubDetId::SubDetEnum det, AnaTrackB* selTracks[]) {
//**************************************************

  int count = 0;

  for (unsigned int it = 0; it < bunch.Particles.size(); it++) {
    AnaTrackB* track = static_cast<AnaTrackB*>(bunch.Particles[it]);
    if (anaUtils::TrackUsesOnlyDet(*track, det)) {
      selTracks[count] = track;
      ++count;
    }
  }

  return count;
}


//**************************************************
int anaUtils::GetAllTracksUsingSubdet1AndNoSubdet2(const AnaBunchB& bunch, AnaTrackB* selTracks[]) {
//**************************************************

  int count = 0;
  for (unsigned int it = 0; it < bunch.Particles.size(); it++) {
    AnaTrackB* track = static_cast<AnaTrackB*>(bunch.Particles[it]);
    if (!anaUtils::TrackUsesDet(*track, SubDetId::kSubdet2) && anaUtils::TrackUsesDet(*track, SubDetId::kSubdet1)) {
      selTracks[count] = track;
      ++count;
    }
  }

  // Sort by decreasing number of hits
  std::sort(&selTracks[0], &selTracks[count], AnaParticleB::CompareNHits);

  return count;
}


//**************************************************
int anaUtils::GetAllTracksUsingSubdet1orSubdet2(const AnaBunchB& bunch, AnaTrackB* selTracks[]) {
//**************************************************

  int count = 0;
  for (unsigned int it = 0; it < bunch.Particles.size(); it++) {
    AnaTrackB* track = static_cast<AnaTrackB*>(bunch.Particles[it]);
    if ( ! (anaUtils::TrackUsesDet(*track, SubDetId::kSubdet2) || anaUtils::TrackUsesDet(*track, SubDetId::kSubdet1))) {
      selTracks[count] = track;
      ++count;
    }
  }

  return count;
}



