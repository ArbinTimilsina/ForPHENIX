


void jetAntiKTcode(){

/////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Make anti-kt jet here
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  float R = 0.3;
  fastjet::JetDefinition *antikt = new fastjet::JetDefinition(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);

  std::vector<fastjet::PseudoJet> jetParticles_all;
  jetParticles_all.clear();

  unsigned int indexTotal = tracks_list.size() + clusters_list.size();
  unsigned int indexCharged = tracks_list.size();

  float particlePt[indexTotal];
  fill(particlePt, particlePt + indexTotal / sizeof(float), -999.9);
  bool particleErtTrigger[indexTotal];
  fill(particleErtTrigger, particleErtTrigger + indexTotal / sizeof(bool), false);

  int index = 0;
  for (unsigned int h = 0; h < tracks_list.size(); h++)
    {
      fastjet::PseudoJet pseudoCharged(tracks_list[h].px,
				       tracks_list[h].py,
				       tracks_list[h].pz,
				       tracks_list[h].mom);
      pseudoCharged.set_user_index(index);
      particlePt[index] = tracks_list[h].pT;
      particleErtTrigger[index] = false;
      jetParticles_all.push_back(pseudoCharged);
      index++;
    }

  for (unsigned int n = 0; n < clusters_list.size(); n++)
    {
      fastjet::PseudoJet pseudoNeutral(clusters_list[n].px,
				       clusters_list[n].py,
				       clusters_list[n].pz,
				       clusters_list[n].energy);
      pseudoNeutral.set_user_index(index);
      particlePt[index] = clusters_list[n].pT;
      particleErtTrigger[index] = clusters_list[n].ertFired;
      jetParticles_all.push_back(pseudoNeutral);
      index++;
    }

  fastjet::ClusterSequence jetAll(jetParticles_all, *antikt);
  std::vector<fastjet::PseudoJet> fastAll = jetAll.inclusive_jets();
  for (unsigned int n = 0; n < fastAll.size(); n++)
    {
      fastjet::PseudoJet aFastJet = fastAll[n];

      float chargedPt    = 0.0;

      float jetPt  = aFastJet.perp();
      float jetEta = aFastJet.pseudorapidity();
      float jetPhi = phiReduce(aFastJet.phi());

      vector<fastjet::PseudoJet> constituents = jetAll.constituents(aFastJet);
      unsigned int nconst = constituents.size();

      bool ertTriggerFired = false;
      for (unsigned int iconst = 0; iconst < nconst; iconst++)
	{
	  unsigned int indx = constituents[iconst].user_index();

	  if (indx < indexCharged)//Charged particles
	    {
	      chargedPt += particlePt[indx];
	    }
	  else  //Neutral particles
	    {
	      if(particleErtTrigger[indx] == true)
		{
		  ertTriggerFired = true;
		}
	    }
	}

      int jetArm = 1;
      if (jetPhi > 1.57)
	{
	  jetArm = 0;
	}

      if(!ertTriggerFired)
	{
	  continue;
	}

      float jetCf = chargedPt / jetPt;
      if((jetPt > 6.0) && ((float)nconst >= 3.0))
	{
	  hCf->Fill(jetCf);
	}

      //Jet selection
      bool passJetLevelCuts = (jetPt > 6.0) && ((float)nconst >= 3.0);
      if(!passJetLevelCuts)
	{
	  continue;
	}
      hJets->Fill(jetPt);

      jets temp;
      temp.arm           = jetArm;
      temp.pT            = jetPt;
      temp.eta           = jetEta;
      temp.phi           = jetPhi;
      temp.nc            = (float)nconst;
      temp.cf            = jetCf;
      jets_list.push_back(temp);

      nJets++;
    }
  delete antikt;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////



}//end jetAntiKTcode()
