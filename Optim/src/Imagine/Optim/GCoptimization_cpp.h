
// will leave this one just for the laughs :)
//#define olga_assert(expr) assert(!(expr))

/////////////////////////////////////////////////////////////////////////////////////////////////
//   First we have functions for the base class
/////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor for base class                                                       
template <typename V,typename T>
GCoptimization<V,T>::GCoptimization(SiteID nSites, LabelID nLabels) 
: m_random_label_order(false)
, m_datacostIndividual(0)
, m_smoothcostIndividual(0)
, m_datacostFn(0)
, m_smoothcostFn(0)
, m_giveDataEnergyInternal(0)
, m_giveSmoothEnergyInternal(0)
, m_set_up_n_links_expansion(0)
, m_set_up_t_links_expansion(0)
, m_set_up_n_links_swap(0)
, m_set_up_t_links_swap(0)
, m_datacostFnDelete(0)
, m_smoothcostFnDelete(0)
{
	assert( nLabels > 1 && nSites > 0);
	m_num_labels = nLabels;
	m_num_sites  = nSites;


	m_lookupSiteVar = (SiteID *)  new SiteID[m_num_sites];
	m_labelTable    = (LabelID *) new LabelID[m_num_labels];
	m_labeling      = (LabelID *) new LabelID[m_num_sites];

	if ( !m_lookupSiteVar || !m_labelTable || !m_labeling){
		if (m_lookupSiteVar) delete [] m_lookupSiteVar;
		if (m_labelTable) delete [] m_labelTable;
		if (m_labeling) delete [] m_labeling;
		handleError("Not enough memory");

		
	}
	
	for ( LabelID i = 0; i < m_num_labels; i++ )
		m_labelTable[i] = i;
		

	for ( SiteID i = 0; i < m_num_sites; i++ ){
		m_labeling[i] = (LabelID) 0;
		m_lookupSiteVar[i] = (SiteID) -1;
	}
	
	//srand((unsigned int) time(NULL));
}

//-------------------------------------------------------------------

template <typename V,typename T>
GCoptimization<V,T>::~GCoptimization()
{
	delete [] m_labelTable;
	delete [] m_lookupSiteVar;
	delete [] m_labeling;

	if (m_datacostFnDelete) m_datacostFnDelete(m_datacostFn);
	if (m_smoothcostFnDelete) m_smoothcostFnDelete(m_smoothcostFn);

	if (m_datacostIndividual) delete [] m_datacostIndividual;
	if (m_smoothcostIndividual) delete [] m_smoothcostIndividual;
}


//------------------------------------------------------------------
template <typename V,typename T>
void GCoptimization<V,T>::setDataCost(DataCostFn fn) { 
	specializeDataCostFunctor(DataCostFnFromFunction(fn));
}


//------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::setDataCost(DataCostFnExtra fn, void *extraData) { 
	specializeDataCostFunctor(DataCostFnFromFunctionExtra(fn, extraData));
}

//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::setDataCost(EnergyTermType *dataArray) {
	specializeDataCostFunctor(DataCostFnFromArray(dataArray, m_num_labels));
}

//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::setDataCost(SiteID s, LabelID l, EnergyTermType e) {
	if (!m_datacostIndividual) {
		if ( m_datacostFn ) handleError("Data Costs are already initialized");
		m_datacostIndividual = new EnergyTermType[m_num_sites*m_num_labels];
		memset(m_datacostIndividual, 0, m_num_sites*m_num_labels*sizeof(EnergyTermType));
		specializeDataCostFunctor(DataCostFnFromArray(m_datacostIndividual, m_num_labels));
	}
	m_datacostIndividual[s*m_num_labels + l] = e;
}

//-------------------------------------------------------------------
template <typename V,typename T>
void GCoptimization<V,T>::setDataCostFunctor(DataCostFunctor* f) {
	if ( m_datacostFn ) handleError("Data Costs are already initialized");

	m_datacostFn = f;
	m_giveDataEnergyInternal   = &GCoptimization<V,T>::template giveDataEnergyInternal<DataCostFunctor>;
	m_set_up_t_links_expansion = &GCoptimization<V,T>::template set_up_t_links_expansion<DataCostFunctor>;
	m_set_up_t_links_swap      = &GCoptimization<V,T>::template set_up_t_links_swap<DataCostFunctor>;
}


//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::setSmoothCost(SmoothCostFn fn) {
	specializeSmoothCostFunctor(SmoothCostFnFromFunction(fn));
}

//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::setSmoothCost(SmoothCostFnExtra fn, void* extraData) {
	specializeSmoothCostFunctor(SmoothCostFnFromFunctionExtra(fn, extraData));
}

//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::setSmoothCost(EnergyTermType *smoothArray) {
	specializeSmoothCostFunctor(SmoothCostFnFromArray(smoothArray, m_num_labels));
}

//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::setSmoothCost(LabelID l1, LabelID l2, EnergyTermType e){
	if (!m_smoothcostIndividual) {
		if ( m_smoothcostFn ) handleError("Smoothness Costs are already initialized");
		m_smoothcostIndividual = new EnergyTermType[m_num_labels*m_num_labels];
		memset(m_smoothcostIndividual, 0, m_num_labels*m_num_labels*sizeof(EnergyTermType));
		specializeSmoothCostFunctor(SmoothCostFnFromArray(m_smoothcostIndividual, m_num_labels));
	} 
	m_smoothcostIndividual[l1*m_num_labels + l2] = e;
}

//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::setSmoothCostFunctor(SmoothCostFunctor* f) {
	if ( m_smoothcostFn ) handleError("Smoothness Costs are already initialized");

	m_smoothcostFn = f;
	m_giveSmoothEnergyInternal = &GCoptimization<V,T>::template giveSmoothEnergyInternal<SmoothCostFunctor>;
	m_set_up_n_links_expansion = &GCoptimization<V,T>::template set_up_n_links_expansion<SmoothCostFunctor>;
	m_set_up_n_links_swap      = &GCoptimization<V,T>::template set_up_n_links_swap<SmoothCostFunctor>;
}

//-------------------------------------------------------------------
 
template <typename V,typename T>
typename GCoptimization<V,T>::EnergyType GCoptimization<V,T>::giveSmoothEnergy()
{
	return( (this->*m_giveSmoothEnergyInternal)());
}
//-------------------------------------------------------------------
 
template <typename V,typename T>
typename GCoptimization<V,T>::EnergyType GCoptimization<V,T>::giveDataEnergy()
{
	return( (this->*m_giveDataEnergyInternal)());
}

//-------------------------------------------------------------------
 
template <typename V,typename T>
typename GCoptimization<V,T>::EnergyType GCoptimization<V,T>::computeEnergy()
{
	if (readyToOptimise())
		return( (this->*m_giveDataEnergyInternal)()+ (this->*m_giveSmoothEnergyInternal)());
	else{ 
		handleError("Not ready to optimize yet. Set up data and smooth costs first");
		return(0);
	}
}

//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::scramble_label_table()
{
   LabelID r1,r2,temp;
   LabelID num_times,cnt;

   num_times = m_num_labels;

   for ( cnt = 0; cnt < num_times; cnt++ )
   {
      r1 = rand()%m_num_labels;  
      r2 = rand()%m_num_labels;  

      temp             = m_labelTable[r1];
      m_labelTable[r1] = m_labelTable[r2];
      m_labelTable[r2] = temp;
   }
}


//-------------------------------------------------------------------

template <typename V,typename T>
typename GCoptimization<V,T>::EnergyType GCoptimization<V,T>::expansion(int max_num_iterations )
{
	int curr_cycle = 1;
	EnergyType new_energy,old_energy;
	

	new_energy = computeEnergy();

	old_energy = new_energy+1;

	while ( old_energy > new_energy  && curr_cycle <= max_num_iterations)
	{
		old_energy = new_energy;
		new_energy = oneExpansionIteration();
		curr_cycle++;	
	}

	return(new_energy);
}

//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::setLabelOrder(bool RANDOM_LABEL_ORDER)
{
	m_random_label_order = RANDOM_LABEL_ORDER;
}

//-------------------------------------------------------------------

template <typename V,typename T>
bool GCoptimization<V,T>::readyToOptimise()
{
	if (!m_smoothcostFn)
	{
		handleError("Smoothness term is not set up yet!");
		return(false);
	}
	if (!m_datacostFn)
	{
		handleError("Data term is not set up yet!");
		return(false);
	}
	return(true);

}

//------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::handleError(const char *message)
{
	throw  GCException(message);
}


//-------------------------------------------------------------------//
//                  METHODS for EXPANSION MOVES                      //  
//-------------------------------------------------------------------//

template <typename V,typename T>
void GCoptimization<V,T>::set_up_expansion_energy(SiteID size, LabelID alpha_label,Energy<V,T> *e,
											VarID *variables, SiteID *activeSites )
{

	
	(this->*m_set_up_t_links_expansion)(size,alpha_label,e,variables,activeSites);
	(this->*m_set_up_n_links_expansion)(size,alpha_label,e,variables,activeSites);
}


//-------------------------------------------------------------------
// Sets up the energy and optimizes it. Also updates labels of sites, if needed. ActiveSites
// are the sites participating in expansion. They are not labeled alpha currrently                                                                    
template <typename V,typename T>
void GCoptimization<V,T>::solveExpansion(SiteID size,SiteID *activeSites,LabelID alpha_label)
{
	SiteID i,site;
	
	if ( !readyToOptimise() ) handleError("Set up data and smoothness terms first. ");
	if ( size == 0 ) return;

	
	Energy<V,T> *e = new Energy<V,T>();

	VarID *variables = (VarID *) new VarID[size];

	for ( i = 0; i < size; i++ )
		variables[i] = e ->add_variable();

	set_up_expansion_energy(size,alpha_label,e,variables,activeSites);
		
	/*typename Energy<V,T>::TotalValue Emin =*/ e -> minimize();
		
	for ( i = 0; i < size; i++ )
	{
		site = activeSites[i];
		if ( m_labeling[site] != alpha_label )
		{
			if ( e->get_var(variables[i]) == 0 )
			{
				m_labeling[site] = alpha_label;
			}
		}
		m_lookupSiteVar[site] = -1;
	}

	delete [] variables;
	delete e;
}
//-------------------------------------------------------------------
// alpha expansion on all sites not currently labeled alpha                         
template <typename V,typename T>
void GCoptimization<V,T>::alpha_expansion(LabelID alpha_label)
{
	SiteID i  = 0, size = 0;
	assert( alpha_label >= 0 && alpha_label < m_num_labels);

	SiteID *activeSites = new SiteID[m_num_sites];
	
	for ( i = 0; i < m_num_sites; i++ )
	{
		if ( m_labeling[i] != alpha_label )
		{
			activeSites[size] = i;
			m_lookupSiteVar[i] = size;
			size++;
		}
	}

	solveExpansion(size,activeSites,alpha_label);

	delete [] activeSites;
}

//-------------------------------------------------------------------
// alpha expansion on subset of sites in array *sites                             
template <typename V,typename T>
void GCoptimization<V,T>::alpha_expansion(LabelID alpha_label, SiteID *sites, SiteID num )
{
	SiteID i,size = 0; 
	SiteID *activeSites = new SiteID[num];
	
	for ( i = 0; i < num; i++ )
	{
		if ( m_labeling[sites[i]] != alpha_label )
		{
			m_lookupSiteVar[sites[i]] = size;
			activeSites[size] = sites[i];
			size++;
		}
	}

	solveExpansion(size,activeSites,alpha_label);
	
	delete [] activeSites;
}

//-------------------------------------------------------------------

template <typename V,typename T>
typename GCoptimization<V,T>::EnergyType GCoptimization<V,T>::oneExpansionIteration()
{
	LabelID next;

	if (m_random_label_order) scramble_label_table();
	
	for (next = 0;  next < m_num_labels;  next++ )
	{
		alpha_expansion(m_labelTable[next]);
	}
	
	return(computeEnergy());
}

//-------------------------------------------------------------------//
//                  METHODS for SWAP MOVES                           //  
//-------------------------------------------------------------------//

template <typename V,typename T>
typename GCoptimization<V,T>::EnergyType GCoptimization<V,T>::swap(int max_num_iterations )
{
	
	int curr_cycle = 1;
	EnergyType new_energy,old_energy;
	

	new_energy = computeEnergy();

	old_energy = new_energy+1;

	while ( old_energy > new_energy  && curr_cycle <= max_num_iterations)
	{
		old_energy = new_energy;
		new_energy = oneSwapIteration();
		
		curr_cycle++;	
	}

	return(new_energy);
}

//--------------------------------------------------------------------------------

template <typename V,typename T>
typename GCoptimization<V,T>::EnergyType GCoptimization<V,T>::oneSwapIteration()
{
	LabelID next,next1;
   
	if (m_random_label_order) scramble_label_table();
		

	for (next = 0;  next < m_num_labels;  next++ )
		for (next1 = m_num_labels - 1;  next1 >= 0;  next1-- )
			if ( m_labelTable[next] < m_labelTable[next1] )
				alpha_beta_swap(m_labelTable[next],m_labelTable[next1]);

	
	return(computeEnergy());
}

//---------------------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::alpha_beta_swap(LabelID alpha_label, LabelID beta_label)
{
	assert( alpha_label >= 0 && alpha_label < m_num_labels && beta_label >= 0 && beta_label < m_num_labels);

	SiteID i  = 0, size = 0;
	SiteID *activeSites = new SiteID[m_num_sites];
	
	for ( i = 0; i < m_num_sites; i++ )	{
		if ( m_labeling[i] == alpha_label || m_labeling[i] == beta_label ){
			activeSites[size] = i;
			m_lookupSiteVar[i] = size;
			size++;
		}
	}

	solveSwap(size,activeSites,alpha_label,beta_label);

	delete [] activeSites;

}
//-----------------------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::alpha_beta_swap(LabelID alpha_label, LabelID beta_label, 
		                   SiteID *alphaSites, SiteID alpha_size, SiteID *betaSites, 
						   SiteID beta_size)

{
	assert( !(alpha_label < 0 || alpha_label >= m_num_labels || beta_label < 0 || beta_label >= m_num_labels) );
	SiteID i,site,size = 0;
		
	SiteID *activeSites = new SiteID[alpha_size+beta_size];

	for ( i = 0; i < alpha_size; i++ )
	{
		site = alphaSites[i];

		assert(site >= 0 && site < m_num_sites );
		assert(m_labeling[site] == alpha_label);

		activeSites[size]     = site;
		m_lookupSiteVar[site] = size;
		size++;

	}

	for ( i = 0; i < beta_size; i++ )
	{
		site = betaSites[i];
		assert(site >= 0 && site < m_num_sites );
		assert(m_labeling[site] == beta_label);

		activeSites[size]     = site;
		m_lookupSiteVar[site] = size;
		size++;

	}

	solveSwap(size,activeSites,alpha_label,beta_label);

	delete [] activeSites;
}

//-----------------------------------------------------------------------------------
template <typename V,typename T>
void GCoptimization<V,T>::solveSwap(SiteID size,SiteID *activeSites,LabelID alpha_label,
							   LabelID beta_label)
{

	SiteID i,site;
	
	if ( !readyToOptimise() ) handleError("Set up data and smoothness terms first. ");
	if ( size == 0 ) return;

	
	Energy<V,T> *e = new Energy<V,T>();
	VarID *variables = (VarID *) new VarID[size];

	for ( i = 0; i < size; i++ )
		variables[i] = e ->add_variable();

	set_up_swap_energy(size,alpha_label,beta_label,e,variables,activeSites);
		
	typename Energy<V,T>::TotalValue Emin = e -> minimize();
		
	for ( i = 0; i < size; i++ )
	{
		site = activeSites[i];
		if ( e->get_var(variables[i]) == 0 )
			m_labeling[site] = alpha_label;
		else m_labeling[site] = beta_label;
		m_lookupSiteVar[site] = -1;
	}

	delete [] variables;
	delete e;

}
//-----------------------------------------------------------------------------------

template <typename V,typename T>
void GCoptimization<V,T>::set_up_swap_energy(SiteID size,LabelID alpha_label,LabelID beta_label,
										Energy<V,T>* e,typename Energy<V,T>::Var *variables,SiteID *activeSites)
{
    (this->*m_set_up_t_links_swap)(size,alpha_label,beta_label,e,variables,activeSites);
	(this->*m_set_up_n_links_swap)(size,alpha_label,beta_label,e,variables,activeSites);
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for the GCoptimizationGridGraph, derived from GCoptimization
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename V,typename T>
GCoptimizationGridGraph<V,T>::GCoptimizationGridGraph(SiteID width, SiteID height,LabelID num_labels)
						:GCoptimization<V,T>(width*height,num_labels)
{
	
	assert( (width > 1) && (height > 1) && (num_labels > 1 ));

	m_weightedGraph = 0;

	for (int  i = 0; i < 4; i ++ )	m_unityWeights[i] = 1;

	m_width  = width;
	m_height = height;

	m_numNeighbors = new SiteID[this->m_num_sites];
	m_neighbors = new SiteID[4*this->m_num_sites];

	SiteID indexes[4] = {-1,1,-m_width,m_width};

	SiteID indexesL[3] = {1,-m_width,m_width};
	SiteID indexesR[3] = {-1,-m_width,m_width};
	SiteID indexesU[3] = {1,-1,m_width};
	SiteID indexesD[3] = {1,-1,-m_width};

	SiteID indexesUL[2] = {1,m_width};
	SiteID indexesUR[2] = {-1,m_width};
	SiteID indexesDL[2] = {1,-m_width};
	SiteID indexesDR[2] = {-1,-m_width};
    
	setupNeighbData(1,m_height-1,1,m_width-1,4,indexes);

	setupNeighbData(1,m_height-1,0,1,3,indexesL);
	setupNeighbData(1,m_height-1,m_width-1,m_width,3,indexesR);
	setupNeighbData(0,1,1,width-1,3,indexesU);
	setupNeighbData(m_height-1,m_height,1,m_width-1,3,indexesD);

	setupNeighbData(0,1,0,1,2,indexesUL);
	setupNeighbData(0,1,m_width-1,m_width,2,indexesUR);
	setupNeighbData(m_height-1,m_height,0,1,2,indexesDL);
	setupNeighbData(m_height-1,m_height,m_width-1,m_width,2,indexesDR);
}

//-------------------------------------------------------------------

template <typename V,typename T>
GCoptimizationGridGraph<V,T>::~GCoptimizationGridGraph()
{
	delete [] m_numNeighbors;
	delete [] m_neighbors;
	if (m_weightedGraph) delete [] m_neighborsWeights;
}


//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimizationGridGraph<V,T>::setupNeighbData(SiteID startY,SiteID endY,SiteID startX,
											  SiteID endX,SiteID maxInd,SiteID *indexes)
{
	SiteID x,y,pix;
	SiteID n;

	for ( y = startY; y < endY; y++ )
		for ( x = startX; x < endX; x++ )
		{
			pix = x+y*m_width;
			m_numNeighbors[pix] = maxInd;

			for (n = 0; n < maxInd; n++ )
				m_neighbors[pix*4+n] = pix+indexes[n];
		}
}

//-------------------------------------------------------------------

template <typename V,typename T>
bool GCoptimizationGridGraph<V,T>::readyToOptimise()
{
	GCoptimization<V,T>::readyToOptimise();
	return(true);
}

//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimizationGridGraph<V,T>::setSmoothCostVH(EnergyTermType *smoothArray, EnergyTermType *vCosts, EnergyTermType *hCosts)
{
	setSmoothCost(smoothArray);
	m_weightedGraph = 1;
	computeNeighborWeights(vCosts,hCosts);
}

//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimizationGridGraph<V,T>::giveNeighborInfo(SiteID site, SiteID *numSites, SiteID **neighbors, EnergyTermType **weights)
{
	*numSites  = m_numNeighbors[site];
	*neighbors = &m_neighbors[site*4];
	
	if (m_weightedGraph) *weights  = &m_neighborsWeights[site*4];
	else *weights = m_unityWeights;
}

//-------------------------------------------------------------------

template <typename V,typename T>
void GCoptimizationGridGraph<V,T>::computeNeighborWeights(EnergyTermType *vCosts,EnergyTermType *hCosts)
{
	SiteID i,n,nSite;
	typename GCoptimization<V,T>::EnergyTermType weight;

	
	m_neighborsWeights = new EnergyTermType[this->m_num_sites*4];

	for ( i = 0; i < this->m_num_sites; i++ )
	{
		for ( n = 0; n < m_numNeighbors[i]; n++ )
		{
			nSite = m_neighbors[4*i+n];
			if ( i-nSite == 1 )            weight = hCosts[nSite];
			else if (i-nSite == -1 )       weight = hCosts[i];
			else if ( i-nSite == m_width ) weight = vCosts[nSite];
			else if (i-nSite == -m_width ) weight = vCosts[i];
	
			m_neighborsWeights[i*4+n] = weight;
		}
	}

}
////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for the GCoptimizationGeneralGraph, derived from GCoptimization
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename V,typename T>
GCoptimizationGeneralGraph<V,T>::GCoptimizationGeneralGraph(SiteID num_sites,LabelID num_labels):GCoptimization<V,T>(num_sites,num_labels)
{
	assert( num_sites > 1 && num_labels > 1 );
    
	m_neighborsIndexes = 0;
	m_neighborsWeights = 0;
	m_numNeighbors     = 0;
	m_neighbors        = 0;

	m_needTodeleteNeighbors        = true;
	m_needToFinishSettingNeighbors = true;
}

//------------------------------------------------------------------

template <typename V,typename T>
GCoptimizationGeneralGraph<V,T>::~GCoptimizationGeneralGraph()
{
	if (m_needToFinishSettingNeighbors) delete [] m_neighbors;

	if ( m_numNeighbors && m_needTodeleteNeighbors )
	{
		for ( SiteID i = 0; i < this->m_num_sites; i++ )
		{
			if (m_numNeighbors[i] != 0 )
			{
				delete [] m_neighborsIndexes[i];
				delete [] m_neighborsWeights[i];
			}
		}

		delete [] m_numNeighbors;
		delete [] m_neighborsIndexes;
		delete [] m_neighborsWeights;
	}

	
}

//-------------------------------------------------------------------

template <typename V,typename T>
bool GCoptimizationGeneralGraph<V,T>::readyToOptimise()
{
	GCoptimization<V,T>::readyToOptimise();
	if (m_needToFinishSettingNeighbors) finishSettingNeighbors();
	return(true);
}

//------------------------------------------------------------------

template <typename V,typename T>
void GCoptimizationGeneralGraph<V,T>::finishSettingNeighbors()
{
	
	Neighbor *tmp;
	SiteID i,site,count;

	m_needToFinishSettingNeighbors = false;
	
	EnergyTermType *tempWeights = new EnergyTermType[this->m_num_sites];
	SiteID *tempIndexes       = new SiteID[this->m_num_sites];
	
	if ( !tempWeights || !tempIndexes ) this->handleError("Not enough memory");

	m_numNeighbors     = new SiteID[this->m_num_sites];
	m_neighborsIndexes = new SiteID*[this->m_num_sites];
	m_neighborsWeights = new EnergyTermType*[this->m_num_sites];
	
	if ( !m_numNeighbors || !m_neighborsIndexes || !m_neighborsWeights ) this->handleError("Not enough memory");

	for ( site = 0; site < this->m_num_sites; site++ )
	{
		if ( !m_neighbors[site].isEmpty() )
		{
			m_neighbors[site].setCursorFront();
			count = 0;
			
			while ( m_neighbors[site].hasNext() )
			{
				tmp = (Neighbor *) (m_neighbors[site].next());
				tempIndexes[count] =  tmp->to_node;
				tempWeights[count] =  tmp->weight;
				count++;
			}
			m_numNeighbors[site]     = count;
			m_neighborsIndexes[site] = new SiteID[count];
			m_neighborsWeights[site] = new EnergyTermType[count];
			
			if ( !m_neighborsIndexes[site] || !m_neighborsWeights[site] ) this->handleError("Not enough memory");
			
			for ( i = 0; i < count; i++ )
			{
				m_neighborsIndexes[site][i] = tempIndexes[i];
				m_neighborsWeights[site][i] = tempWeights[i];
			}
		}
		else m_numNeighbors[site] = 0;

	}

	delete [] tempIndexes;
	delete [] tempWeights;
	delete [] m_neighbors;
	

}
//------------------------------------------------------------------------------

template <typename V,typename T>
void GCoptimizationGeneralGraph<V,T>::giveNeighborInfo(SiteID site, SiteID *numSites, 
												  SiteID **neighbors, EnergyTermType **weights)
{
	(*numSites)  =  m_numNeighbors[site];
	(*neighbors) = m_neighborsIndexes[site];
	(*weights)   = m_neighborsWeights[site];
}


//------------------------------------------------------------------

template <typename V,typename T>
void GCoptimizationGeneralGraph<V,T>::setNeighbors(SiteID site1, SiteID site2, EnergyTermType weight)
{

	assert( site1 < this->m_num_sites && site1 >= 0 && site2 < this->m_num_sites && site2 >= 0);
	if ( m_needToFinishSettingNeighbors == false ) this->handleError("Already set up neighborhood system");

	if ( !m_neighbors )
	{
		m_neighbors = (LinkedBlockList *) new LinkedBlockList[this->m_num_sites];
		if ( !m_neighbors ) this->handleError("Not enough memory");
	}

	Neighbor *temp1 = (Neighbor *) new Neighbor;
	Neighbor *temp2 = (Neighbor *) new Neighbor;

	temp1->weight  = weight;
	temp1->to_node = site2;

	temp2->weight  = weight;
	temp2->to_node = site1;

	m_neighbors[site1].addFront(temp1);
	m_neighbors[site2].addFront(temp2);
	
}
//------------------------------------------------------------------

template <typename V,typename T>
void GCoptimizationGeneralGraph<V,T>::setAllNeighbors(SiteID *numNeighbors,SiteID **neighborsIndexes,
												 EnergyTermType **neighborsWeights)
{
	m_needTodeleteNeighbors = false;
	m_needToFinishSettingNeighbors = false;
	if ( m_neighbors ) this->handleError("Already set up neighborhood system");

	m_numNeighbors     = numNeighbors;
	m_neighborsIndexes = neighborsIndexes;
	m_neighborsWeights = neighborsWeights;
}
