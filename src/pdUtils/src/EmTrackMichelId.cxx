#include "EmTrackMichelId.hxx"

/*nnet::EmTrackMichelId::EmTrackMichelId(EmTrackMichelId::Parameters const& config) :
        EDProducer{config},
        fBatchSize(config().BatchSize()),
        fPointIdAlg(config().PointIdAlg()),
        fMVAWriter(producesCollector(), "emtrkmichel"),
        fWireProducerLabel(config().WireLabel()),
        fHitModuleLabel(config().HitModuleLabel()),
        fClusterModuleLabel(config().ClusterModuleLabel()),
        fTrackModuleLabel(config().TrackModuleLabel()),
        fViews(config().Views()),

        fNewClustersTag(
            config.get_PSet().get<std::string>("module_label"), "",
            art::ServiceHandle<art::TriggerNamesService const>()->getProcessName())
{
    fMVAWriter.produces_using< my_recob::Hit >();

    if (!fClusterModuleLabel.label().empty())
    {
        produces< std::vector<my_recob::Cluster> >();
            produces< art::Assns<my_recob::Cluster, my_recob::Hit> >();

        fMVAWriter.produces_using< my_recob::Cluster >();
        fDoClusters = true;
    }
    else { fDoClusters = false; }

    if (!fTrackModuleLabel.label().empty())
    {
        fMVAWriter.produces_using< my_recob::Track >();
        fDoTracks = true;
    }
    else { fDoTracks = false; }
    }*/

nnet::EmTrackMichelId::EmTrackMichelId():
//        EDProducer{config},
        fBatchSize(1),
        //        fPointIdAlg(config().PointIdAlg()),
        //        fMVAWriter(producesCollector(), "emtrkmichel"),
        fWireProducerLabel("wclsdatasp:gaus"),
        fHitModuleLabel("hitpdune"),
        fClusterModuleLabel("pandora"),
        fTrackModuleLabel("dummy")  // anselmo
//        fViews(config().Views())  // anselmo
//        fNewClustersTag()
{


  fViews.clear();  //anselmo
  
  if (!fClusterModuleLabel.empty())fDoClusters = true;
  else fDoClusters = false;

  if (!fTrackModuleLabel.empty())fDoTracks = true;
  else fDoTracks = false;

  //fNewClustersTag;? // input tag for the clusters produced by this module

}
// ------------------------------------------------------

void nnet::EmTrackMichelId::produce(AnaEvent & evt)
{


  std::cout << "next event: " << evt.EventInfo.Run << " / " << evt.EventInfo.Event << std::endl;
  /*
  mf::LogVerbatim("EmTrackMichelId") << "next event: " << evt.run() << " / " << evt.id().event();

  auto wireHandle = evt.getValidHandle< std::vector<my_recob::Wire> >(fWireProducerLabel);
  */

  detinfo::DetectorClocksData clockData;
  detinfo::DetectorPropertiesData detProp;

  
  std::vector<my_recob::Wire> wireHandle;
  unsigned int cryo, tpc, view;
  
  // ******************* get and sort hits ********************

  std::vector<AnaHitPD*> hitListHandle;
  std::vector<AnaHitPD*> hitPtrList;
  //  auto hitListHandle = evt.getValidHandle< std::vector<my_recob::Hit> >(fHitModuleLabel);
  //  art::fill_ptr_vector(hitPtrList, hitListHandle);

  std::cout << "anselmo 0: nparts = " << evt.nParticles << std::endl;  
  for (UInt_t i=0;i<evt.nParticles;i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(evt.Particles[i]);
    for (UInt_t j=0;j<part->Hits[2].size();j++){
      hitListHandle.push_back(&(part->Hits[2][j]));
      hitPtrList.push_back(&(part->Hits[2][j]));
    }
  }
  std::cout << "anselmo 1: nhits = " << hitPtrList.size() << std::endl;

  EmTrackMichelId::cryo_tpc_view_keymap hitMap;
  for (auto const& h : hitPtrList)
    {
      view = h->WireID().Plane;
      if (!isViewSelected(view)) continue;
      
      cryo = h->WireID().Cryostat;
      tpc = h->WireID().TPC;
      
      //      hitMap[cryo][tpc][view].push_back(h.key());
      //      hitMap[cryo][tpc][view].push_back((size_t)h);
      // anselmo
      hitMap[cryo][tpc][view].push_back(h);
    }

    // ********************* classify hits **********************
  //    auto hitID = fMVAWriter.initOutputs<my_recob::Hit>(fHitModuleLabel, hitPtrList.size(), fPointIdAlg.outputLabels());

    std::vector< char > hitInFA(hitPtrList.size(), 0); // tag hits in fid. area as 1, use 0 for hits close to the projectrion edges
    for (auto const & pcryo : hitMap)
    {
        cryo = pcryo.first;
        for (auto const & ptpc : pcryo.second)
        {
            tpc = ptpc.first;
            for (auto const & pview : ptpc.second)
            {
                view = pview.first;
                if (!isViewSelected(view)) continue; // should not happen, hits were selected

                fPointIdAlg.setWireDriftData(clockData, detProp, wireHandle, view, tpc, cryo);

                // (1) do all hits in this plane ------------------------------------------------
                for (size_t idx = 0; idx < pview.second.size(); idx += fBatchSize)
                {
                    std::vector< std::pair<unsigned int, float> > points;
                    std::vector< size_t > keys;
                    for (size_t k = 0; k < fBatchSize; ++k)
                    {

                        if (idx + k >= pview.second.size()) { break; } // careful about the tail
                        // anselmo:  use directly the pointer, not using the key
                        //                        size_t h = pview.second[idx+k]; // h is the Ptr< my_recob::Hit >::key()
                        size_t h = (size_t)(pview.second[idx+k]);
                        //                        const AnaHitPD & hit = *(hitPtrList[h]);
                        const AnaHitPD & hit = *(pview.second[idx+k]);
                        points.emplace_back(hit.WireID().Wire, hit.PeakTime());
                        //                        keys.push_back(h);
                        keys.push_back((size_t)h);
                    }
                    auto batch_out = fPointIdAlg.predictIdVectors(points);
                    if (points.size() != batch_out.size())
                    {
                      //                        throw cet::exception("EmTrackMichelId") << "hits processing failed" << std::endl;
                      std::cout << "hits processing failed" << std::endl;
                    }
                    for (size_t k = 0; k < points.size(); ++k)
                    {
                        size_t h = keys[k];
                        //                        fMVAWriter.setOutput(hitID, h, batch_out[k]);
                        if (fPointIdAlg.isInsideFiducialRegion(points[k].first, points[k].second))
                        { hitInFA[h] = 1; }
                    }
                } // hits done ------------------------------------------------------------------
            }
        }
    }

  /*

  // (2) do clusters when hits are ready in all planes ----------------------------------------
    if (fDoClusters)
    {
        // **************** prepare for new clusters ****************
            auto clusters = std::make_unique< std::vector< my_recob::Cluster > >();
            auto clu2hit = std::make_unique< art::Assns< my_recob::Cluster, my_recob::Hit > >();

        // ************** get and sort input clusters ***************
        auto cluListHandle = evt.getValidHandle< std::vector<my_recob::Cluster> >(fClusterModuleLabel);
            std::vector< art::Ptr<my_recob::Cluster> > cluPtrList;
            art::fill_ptr_vector(cluPtrList, cluListHandle);

        EmTrackMichelId::cryo_tpc_view_keymap cluMap;
            for (auto const& c : cluPtrList)
            {
                view = c->Plane().Plane;
                if (!isViewSelected(view)) continue;

                cryo = c->Plane().Cryostat;
                tpc = c->Plane().TPC;

                cluMap[cryo][tpc][view].push_back(c.key());
            }

        auto cluID = fMVAWriter.initOutputs<my_recob::Cluster>(fNewClustersTag, fPointIdAlg.outputLabels());

        unsigned int cidx = 0; // new clusters index
        art::FindManyP< my_recob::Hit > hitsFromClusters(cluListHandle, evt, fClusterModuleLabel);
        std::vector< bool > hitUsed(hitPtrList.size(), false); // tag hits used in clusters
        for (auto const & pcryo : cluMap)
        {
            cryo = pcryo.first;
            for (auto const & ptpc : pcryo.second)
            {
                tpc = ptpc.first;
                for (auto const & pview : ptpc.second)
                {
                    view = pview.first;
                    if (!isViewSelected(view)) continue; // should not happen, clusters were pre-selected

                    for (size_t c : pview.second) // c is the Ptr< my_recob::Cluster >::key()
                    {
                                auto v = hitsFromClusters.at(c);
                                if (v.empty()) continue;

                        for (auto const & hit : v)
                        {
                            if (hitUsed[hit.key()]) { mf::LogWarning("EmTrackMichelId") << "hit already used in another cluster"; }
                            hitUsed[hit.key()] = true;
                        }

                        auto vout = fMVAWriter.getOutput<my_recob::Hit>(v,
                            [&](art::Ptr<my_recob::Hit> const & ptr) { return (float)hitInFA[ptr.key()]; });

                            float pvalue = vout[0] / (vout[0] + vout[1]);
                            mf::LogVerbatim("EmTrackMichelId") << "cluster in tpc:" << tpc << " view:" << view
                            << " size:" << v.size() << " p:" << pvalue;

                                clusters->emplace_back(
                                    my_recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
                                            v.size(), 0.0F, 0.0F, cidx, (geo::View_t)view, v.front()->WireID().planeID()));
                            util::CreateAssn(*this, evt, *clusters, v, *clu2hit);
                            cidx++;

                            fMVAWriter.addOutput(cluID, vout); // add copy of the input cluster
                    }

                    // (2b) make single-hit clusters --------------------------------------------
                    for (size_t h : hitMap[cryo][tpc][view]) // h is the Ptr< my_recob::Hit >::key()
                    {
                        if (hitUsed[h]) continue;

                        auto vout = fMVAWriter.getOutput<my_recob::Hit>(h);
                                float pvalue = vout[0] / (vout[0] + vout[1]);

                                mf::LogVerbatim("EmTrackMichelId") << "single hit in tpc:" << tpc << " view:" << view
                                        << " wire:" << hitPtrList[h]->WireID().Wire << " drift:" << hitPtrList[h]->PeakTime() << " p:" << pvalue;

                                art::PtrVector< my_recob::Hit > cluster_hits;
                                cluster_hits.push_back(hitPtrList[h]);
                                clusters->emplace_back(
                                        my_recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
                                            1, 0.0F, 0.0F, cidx, (geo::View_t)view, hitPtrList[h]->WireID().planeID()));
                                util::CreateAssn(*this, evt, *clusters, cluster_hits, *clu2hit);
                                cidx++;

                                fMVAWriter.addOutput(cluID, vout); // add single-hit cluster tagging unclutered hit
                    }
                    mf::LogVerbatim("EmTrackMichelId") << "...produced " << cidx - pview.second.size() << " single-hit clusters.";
                }
            }
        }

        evt.put(std::move(clusters));
            evt.put(std::move(clu2hit));
        } // all clusters done ----------------------------------------------------------------------
  */
    // (3) do tracks when all hits in all cryo/tpc/plane are done -------------------------------
    if (fDoTracks)
    {
      /*
        auto trkListHandle = evt.getValidHandle< std::vector<my_recob::Track> >(fTrackModuleLabel);
        art::FindManyP< my_recob::Hit > hitsFromTracks(trkListHandle, evt, fTrackModuleLabel);
        std::vector< std::vector< art::Ptr<my_recob::Hit> > > trkHitPtrList(trkListHandle->size());
      */
      std::vector<my_recob::Track>  trkListHandle;
      std::vector<std::vector<AnaHitPD> > hitsFromTracks;
      std::vector<std::vector<AnaHitPD> > trkHitPtrList;


      for (UInt_t i=0;i<evt.nParticles;i++){
        AnaParticlePD* part = static_cast<AnaParticlePD*>(evt.Particles[i]);
        my_recob::Track track(*part);
        trkListHandle.push_back(track);
        hitsFromTracks.push_back(part->Hits[2]);        
      }
      trkHitPtrList.resize(trkListHandle.size());
      
        for (size_t t = 0; t < trkListHandle.size(); ++t)
        {
            auto v = hitsFromTracks.at(t);
            size_t nh[3] = { 0, 0, 0 };
            for (auto const & hptr : v) { ++nh[hptr.View()]; }
            size_t best_view = 2; // collection
            if ((nh[0] >= nh[1]) && (nh[0] > 2 * nh[2])) best_view = 0; // ind1
            if ((nh[1] >= nh[0]) && (nh[1] > 2 * nh[2])) best_view = 1; // ind2
            size_t k = 0;
            while (!isViewSelected(best_view))
            {
                best_view = (best_view + 1) % 3;
                //                if (++k > 3) { throw cet::exception("EmTrackMichelId") << "No views selected at all?" << std::endl; }
                if (++k > 3) { std::cout << "No views selected at all?" << std::endl; }
            }
            for (auto const & hptr : v)
            {
                if (hptr.View() == best_view) trkHitPtrList[t].emplace_back(hptr);
            }
        }
        /*
        auto trkID = fMVAWriter.initOutputs<my_recob::Track>(fTrackModuleLabel, trkHitPtrList.size(), fPointIdAlg.outputLabels());
        for (size_t t = 0; t < trkHitPtrList.size(); ++t) // t is the Ptr< my_recob::Track >::key()
        {
            auto vout = fMVAWriter.getOutput<my_recob::Hit>(trkHitPtrList[t],
                [&](art::Ptr<my_recob::Hit> const & ptr) { return (float)hitInFA[ptr.key()]; });
            fMVAWriter.setOutput(trkID, t, vout);
        }
        */
    }
    // tracks done ------------------------------------------------------------------------------

    //        fMVAWriter.saveOutputs(evt);
}
// ------------------------------------------------------

bool nnet::EmTrackMichelId::isViewSelected(int view) const
{
  if (fViews.empty()) return true;
  else
    {
      bool selected = false;
      for (auto k : fViews) if (k == view) { selected = true; break; }
      return selected;
    }
}
// ------------------------------------------------------