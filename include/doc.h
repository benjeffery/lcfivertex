/*! \page LCIO LCIO Interface
	Description of this package's use of LCIO to store results of processors.
	\section Vertex Storage of Vertexing Result
		The processors ZVTOPZVKINProcessor and ZVTOPZVRESProcessor 
		
		both store their results in the same way. They provide information on the decay vertices found, stored in a collection of LCIO Vertices and a collection of ReconstructedParticles, representing decaying particles that give rise to these vertices, as described in Frank Gaede’s <a href="http://forum.linearcollider.org/index.php?t=tree&goto=548"> Frank Gaede's forum post</a>
		<br>
		In LCIO, Vertices and ReconstructedParticles, are connected as follows:
		<ul>
		<li>Each decay vertex found has a corresponding LCIO::Vertex.</li>
		<li>Each decay vertex found also has a corresponding LCIOReconstructedParticle which respresents the decaying particle and holds kinematic information.</li>
		<li>This accompanying decaying ReconstructedParticle is accessed though Vertex::getAssociatedParticle()</li>
		<li>The descendant tracks which are produced by the particle are attached to the decaying ReconstructedParticle and accessed through ReconstructedParticle::getParticles()</li>
		<li>Each ReconstructedParticle points to its start and end vertex (if any) through ReconstructedParticle::getStartVertex() and ReconstructedParticle::getEndVertex()</li>
		</ul>
		In essence we end up with three types of objects with links between them; Vertices, Decaying ReconstructedParticles, Stable ReconstructedParticles.
		<br>
		Note that the information stored in the ReconstructedParticles created using ZVTOP information differs from the one of those created by the track reconstruction code. The getStartVertex and getEndVertex methods permit building up a representation of a heavy flavour decay chain. In the current ZVTOP version, this is reconstructed as follows: vertices are sorted by increasing radius from the IP, the start vertex of the accompanying ReconstructedParticle is set to point to the previous vertex in the collection in this order. Note that this is only an approximation of the decay chain in an actual physics event, which may differ or be more complex (e.g. one decay vertex may give rise to two instable particles. If correctly reconstructed, their ReconstructedParticles would point to the same start vertex. Of the corresponding ZVTOP vertices, only the one closer to the IP will point to its correct start vertex, while the further one will point to the nearer one instead.
		<br>
		As each RP points to its vertices, to store more than one vertexing result it is necessary to take a copy of the set of RPs that are in each vertex. Each copy points to the original unique RP through ReconstructedParticle::getParticles()[0].
		<br>
		A master ReconstructedParticle object is created that points to all decaying and stable ReconstructedParticles in the decay chain though ReconstructedParticle::getParticles(). The ReconstructedParticle::getStartVertex() is set to the first Vertex in the decay chain (usually the IP). This is the main object for accessing the decay chains as it allows one to iterate over all the ReconstructedParticles contained within.
		<br>
		These objects are then stored in three collections: <ul>
		<li>
		DecayChainRPTracksCollectionName: default name: ZV***DecayChainRPTracks, stores decaying RPs and copies of input RPs for all decay chains.
		</li>
		<li>
		VertexCollection: default name: ZV***Vertices, stores the Vertices for all decay chains.
		</li>
		<li>
		DecayChainCollectionName: default name: ZV***DecayChains, stores the master ReconstructedParticle DecayChainCollectionName, in the same order as the Jets that were input to the algoithm in JetRPCollection
		</li>
		</ul>
		<h4>Example</h4>
		If you wanted to print out the d0 of each LCIO::Track in each Vertex of a decay chain found for the second jet in the jet collection:
		<p>
		\code
		//Get the decay chain master RP collection and get the second decay chain
		LCCollection* DecayChainRPCol = evt->getCollection( _DecayChainRPColName );
		ReconstructedParticle* DecayChainRP = dynamic_cast<ReconstructedParticle*>(DecayChainRPCol->getElementAt(1)
		//Make a list of vertices
		vector<lcio::Vertex*> LCIOVertices;
		//Add the primary first
		LCIOVertices.push_back(DecayChainRP->getStartVertex());
		//Loop over RPs to find all the other vertices
		vector<ReconstructedParticle*> RPs = DecayChainRP->getParticles();
		for (vector<ReconstructedParticle*>::const_iterator iRP = RPs.begin();iRP < RPs.end();++iRP)
		{
			 lcio::Vertex* MyVertex = (*iRP)->getStartVertex();
			 if(MyVertex)
			 {
				 vector<lcio::Vertex*>::const_iterator it = find(LCIOVertices.begin(),LCIOVertices.end(),MyVertex);
				 if (it == LCIOVertices.end())
				 {
					 LCIOVertices.push_back(MyVertex);
				 }
			 }
		}
		//Loop over the vertices
		for (size_t i = 0; i < LCIOVertices.size(); ++i)
		{
			std::cout << "Vertex " << i << "has d0's:" << std::endl;
			
			//Get the Vertex RPs from the Vertices decaying RP
			std::vector<ReconstuctedParticle*> VertexRPs = LCIOVertices[i]->getAssociatedParticle()->getParticles();
			for (size_t j = 0; j < VertexRPs.size(); ++i)
			{
				//Get the Track (remember the Vertex RP is a copy which points to the original)
				std::cout << VertexRPs[j]->getParticles()[0]->getTracks()[0]->getD0() << ",";
			}
			std::cout << std::endl;
		}
		\endcode
		</p>
		<h4>
		Storage of track chi squareds
		</h4>
		If OutputTrackChi2 is set to true the vertexing processors will output the chi squared each track contributes to its vertex.
		<br>
		This is stored in a collection named: TrackRPCollectionName+"TrackChiSquareds" ie. The name of the first collection appended with "TrackChiSquareds"
		<br>
		The chi squareds are stored as a collection of LCFloatVecs in the same order as the tracks in DecayChainRPTracksCollection.
		Currently the LCFloatVec contains one value - the chi squared in the start vertex, but the end vertex could be supported if needed as a second value.
	
		\section FT Storage of Flavour Tag Inputs and Flavour Tag Result
		The lists of variables produced by the Flavour Tag and Inputs are stored in collections of LCFloatVecs, with one LCFloatVec per jet. Within the LCFloatVec, the names and order of the variables are stored in the parameter of the run header as a string vector stored under the name of the LCFloatVec collection. Note that the same LCFloatVec is used for the flavour tag inputs, the vertex charge and further additional variables.
		<h4>Example</h4>
		Get the secondary vertex probability for the second jet:
		<p>
		\code
		//Get the variable names in the FlavourTagInputsCollection
		std::vector<std::string> VarNames;
		(pRun->parameters()).getStringVals(_FlavourTagInputsCollectionName,VarNames);
		//Convert this to a convenient map
		std::map<std::string,unsigned int> IndexOf;
		for (size_t i = 0;i < VarNames.size();++i)
		{
			IndexOf[VarNames[i]] = i;
		}
		//Get the inputs for the second jet
		lcio::LCCollection* pInputs=pEvent->getCollection( _FlavourTagInputsCollectionName );
		LCFloatVec* FTInputs = dynamic_cast<lcio::LCFloatVec*>( pInputs->getElementAt(1) );
		</pre>
		Get the Secondary Probability by indexing the LCFloatVec
		<pre>
		double SecProb = (*FTInputs)[IndexOf["SecondaryVertexProbability"]];
		\endcode
		</p>
		<em>Ben Jeffery - b.jeffery1@physics.ox.ac.uk</em>
*/


/*! \page Cuts Description of Track Selection Cuts
The track selection cuts are mainly from LC note LC-PHSM-2000-021 with minor changes.
<br>These have not yet been optimised for full reconstruction.
	\section IPCuts IP Fitting Cuts
	Details to follow - for now see steering file ipfit.xml
	\section ZVTOPCuts ZVTOP Cuts
	Details to follow - for now see steering file zvres.xml
	\section FlavourTagInputsCuts FlavourTagInputs Cuts
	Details to follow - for now see steering file fti.xml
*/












