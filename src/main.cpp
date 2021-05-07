#include "GateNode.hpp"
#include "QASMparser.h"
#include "util.cpp"
#include <cassert>
#include <cstring>
#include <iostream>
using namespace std;

//adding the below include for unordered_set for closed set implementation and for built in NlogN sort on vectors
//adding unordered_map for multimap implementation of quick lookup of edge neighbors
#include <unordered_set>
//for fast vector sort
#include <algorithm>
#include <unordered_map>
#include <list>
#include <queue>
#include <map>
#include <deque>

#include <climits>
//for priority queue
#include <functional>

//if smart is defined, we will use a distance matrix instead of going through bfs every time
int SMART = 1;
//if latency is 1, we will try and gracefully insert swaps between maximal partitions
int LATENCY= 1;

int numLogical = 0;

vector <tuple<map<int, int>, set<int>, set<GateNode*>> > maximalMapper( map<int, set<int>> architectureEdges, set<GateNode*> startSet);
bool wentTooFar(map<int,set<int>> circuitEdges, map<int,set<int>>architectureEdges,int circuitMaxDegree, int architectureMaxDegree);
map<int,set<int>> rollBack(queue<pair<int,int>>changes, map<int,set<int>>architectureEdges);
map<int,set<int>> addGateEdge(GateNode* chosen, map<int,set<int>>circuitEdges);
 pair<queue<pair<int,int>>,map<int,set<int>>> pruneArchitectureEdges(map<int, set<int>>architectureEdges, unsigned int circuitMinDegree);
 map<int, int> perfectMapper(map<int, set<int>> architectureEdges, set<GateNode*> startSet);
 queue<map<int, int>> getMaximalMappings(set<GateNode*> startGates, map<int, set<int>> architectureEdges);
 vector<vector<int>> getPhysicalDistancesFloydWarshall(map<int, set<int>> architectureEdges);
 int partialMappingStitchingCost(map<int, int> firstMap, map<int, int> secondMap, vector<vector<int>> distanceMatrix, map<int, set<int>> architectureEdges);
 int partialMappingCost(map<int, int> mapping, set<GateNode*> remainingGates, vector<vector<int>> distanceMatrix, map<int, set<int>> architectureEdges);
 pair<map<int, int>, vector<pair<int, int>>> stitchMappings(pair<map<int, int>, set<int>> mapOne, pair<map<int, int>, set<int>> mapTwo, map<int, set<int>> architectureEdges);
void lazySwapCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates, vector<vector<int>> distanceMatrix);
void betterLazySwapCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates, vector<vector<int>> distanceMatrix);

//global distance matrix
 vector<vector<int>> distanceMatrix;

 struct PriorityData {
	 int cost;
	 vector<pair<int, int>> swapsTaken;
	 map<int, int> initialMapping;
	 map<int, int> transformedMapping;

	 PriorityData(int first, vector<pair<int, int>> second, map<int, int> third, map<int, int> fourth) {
		 cost = first;
		 swapsTaken = second;
		 initialMapping = third;
		 transformedMapping = fourth;
	 }
	
	 //reversing this operator because we want a min queue
	 bool operator<(const struct PriorityData& toCompare) const{
		 return cost > toCompare.cost; 
	 }

 };

 //priority queue struct format for optimal swap pathway finder
  struct SecondPriorityData {
	  //cost of transforming from currentMapping to targetMapping 
	 int cost;
	 //keeping track of all the swaps between logical qubits done
		//vector of vectors of sets of pairs. the sets represent all parallel swaps between mappings, and we might have more than 1
	 vector<vector<set<pair<int, int>>>> swapsTaken;
	 //keeping track of the initial mapping for the entire process
	 map<int, int> initialMapping;
	 //this is the next proposed mapping in the next layer of maximal mappings to stitch to
	 map<int, int>  targetMapping;
	//this is our current transformedMapping in whatever maximal mapping layer we are at
	 map<int, int> transformedMapping;
	//keeping track of which gate set we are at 
	 set<GateNode*> firstGates;
	 //keeping track of each section of end gates
	 vector<set<GateNode*>> endingSections;

	 SecondPriorityData(int first, vector<vector<set<pair<int, int>>>> second, map<int, int> third, map<int, int> fourth, map<int,int> fifth, set<GateNode*> gates, vector<set<GateNode*>> endings) {
		 cost = first;
		 swapsTaken = second;
		 initialMapping = third;
		 targetMapping = fourth;
		 transformedMapping = fifth;
		 firstGates = gates;
		 endingSections = endings;

	 }
	
	 //reversing this operator because we want a min queue
	 bool operator<(const struct SecondPriorityData& toCompare) const{
		 return cost > toCompare.cost; 
	 }

 };




int main(int argc, char** argv) {
	char* qasmFileName = NULL;
	char* couplingMapFileName = NULL;
	int latency1 = 1;
	int latency2 = 1;
	int latencySwp = 1;

	//Parse command-line arguments:
	for (int iter = 1; iter < argc; iter++) {
		if (!strcmp(argv[iter], "-latency")) {
			latency1 = atoi(argv[++iter]);
			latency2 = atoi(argv[++iter]);
			latencySwp = atoi(argv[++iter]);
		}
		else if (!qasmFileName) {
			qasmFileName = argv[iter];
		}
		else if (!couplingMapFileName) {
			couplingMapFileName = argv[iter];
		}
		else {
			assert(false);
		}
	}

	//Build dependency graph for the quantum circuit's gates; put dependency graph's roots into a set
	set<GateNode*> firstGates;
	int numLogicalQubits = -1;
	int numGates = -1;
	buildDependencyGraph(qasmFileName, firstGates, numLogicalQubits, numGates);

	//Parse the coupling map; put edges into a set
	int numPhysicalQubits = 0;
	set<pair<int, int> > couplings;
	buildCouplingMap(couplingMapFileName, couplings, numPhysicalQubits);
	assert(numPhysicalQubits >= numLogicalQubits);


	//student code goes here?

	numLogical = numLogicalQubits;
	//debug
	//cout << "PRINTING COUPLING\n";
	for (auto i = couplings.begin();i != couplings.end();i++) {
		//printf("%d connected to %d\n", (*i).first, (*i).second);
	}

	//creating a multimap for quick lookup of architecture edges
	//multimap<int, int> architectureEdges ;
	map<int,set<int>> architectureEdges;

	//easy calculation of node degrees
	map<int, int> architectureNodeDegrees;

	//building edges and degrees
	for (auto i = couplings.begin();i != couplings.end();i++) {
		if (architectureEdges.find(i->first)==architectureEdges.end()){
			set<int> edges;
			edges.insert(i->second);
			(architectureEdges)[i->first]=edges;
		} else{
			(architectureEdges)[i->first].insert(i->second);
		}

		if (architectureEdges.find(i->second)==architectureEdges.end()){
			set<int> edges;
			edges.insert(i->first);
			(architectureEdges)[i->second]=edges;
		} else{
			(architectureEdges)[i->second].insert(i->first);
		}
		
	}
	
	/*
	//cout << "Printing architecture edges!\n";
	for (auto i = architectureEdges.begin();i != architectureEdges.end();i++) {
		////printf("Physical Edge between physical qubit %d and physical qubit %d\n", i->first, i->second);
	}

	//cout << "Printing architecture degrees!\n";
	for (auto i = architectureNodeDegrees.begin();i != architectureNodeDegrees.end();i++) {
	//	//printf("Physical qubit %d has degree %d\n", i->first, i->second);
	}
	*/

	//getting floyd warshall distance matrix
	distanceMatrix = getPhysicalDistancesFloydWarshall(architectureEdges);

	/*TEST FOR DISTANCE CORRECTNESS*/
	/*
	//printf("CONFIRMING DISTANCE MATRIX!\n");
	for (int i = 0;i<numPhysicalQubits;i++) {
		for (int j = 0;j < numPhysicalQubits;j++) {
			//printf("PHYSICAL QUBIT %d is %d distance away from PHYSICAL QUBIT %d\n", i, distanceMatrix[i][j], j);
		}
	}
	*/
			
	/* TEST FOR PERFECT MAPPER */

	/*
	map<int, int> result = perfectMapper(architectureEdges, firstGates);
	if (result.size() == 0) {
		//printf("no perfect mapping found!\n");
	}
	else {
		//printf("location of qubits: ");
		for (int i = 0;i < numLogicalQubits;i++) {
			if (result.find(i) != result.end()) {
				//printf("%d, ", result[i]);
			}
			else {
				//printf("-1, ");
			}
		}
		//printf("\n");
	}
	*/

	/* TEST FOR ALL MAXIMAL MAPPINGS*/
	/*
	queue<map<int, int>> maximalMappings = getMaximalMappings(firstGates, architectureEdges);
	map<int, int> lastMapping;
	map<int, int> initialMapping;
	int count = 0;
	while (maximalMappings.size() > 0) {
		
		map<int, int> mapping = maximalMappings.front();
		maximalMappings.pop();
		if (count == 0) {
			initialMapping = mapping;
		}
		cout << "MAXIMAL MAPPING \n";
		for (int i = 0;i < numLogicalQubits;i++) {
			if (mapping.find(i) != mapping.end()) {
				printf("%d, ", mapping[i]);
			}
			else {
				printf("-1, ");
			}
		}
		printf("\n");
		if (lastMapping.size() == 0) {
			lastMapping = mapping;
		}
		else {
			//getting swaps
			printf("DOING SWAPS!\n");
			set<int> dummy;
			pair<map<int,int>,vector<pair<int, int>>> swapFinderResult = stitchMappings(make_pair(lastMapping,dummy),make_pair(mapping,dummy),architectureEdges);

	
			map<int, int> newInitialMapping = swapFinderResult.first;
			vector<pair<int, int>> swaps = swapFinderResult.second;
			for (auto i = swaps.begin();i != swaps.end();i++) {
				printf("SWAP BETWEEN %d and %d\n", i->first, i->second);
			}
			for (auto z = newInitialMapping.begin();z != newInitialMapping.end();z++) {
				if (initialMapping.find(z->first) == initialMapping.end()) {
					initialMapping[z->first] = z->second;
				}
			}
		}
	}
		
	*/	

	/*Test for swaps*/
	/*
	map<int, int> mapOne;
	set<int> mappedOne;
	mapOne[0] = 2;
	mappedOne.insert(2);
	mapOne[1] = 0;
	mappedOne.insert(0);
	mapOne[2] = 1;
	mappedOne.insert(1);
	mapOne[3] = 3;
	mappedOne.insert(3);
	map<int, int> mapTwo;
	set<int> mappedTwo;
	mapTwo[1] = 0;
	mappedTwo.insert(0);
	mapTwo[2] = 1;
	mappedTwo.insert(0);
	mapTwo[3] = 2;
	mappedTwo.insert(0);
	
	

	pair<map<int,int>,vector<pair<int, int>>> swapFinderResult = stitchMappings(make_pair(mapOne,mappedOne),make_pair(mapTwo,mappedTwo),architectureEdges);
	
	map<int, int> newInitialMapping = swapFinderResult.first;
	vector<pair<int, int>> swaps = swapFinderResult.second;
	for (auto i = swaps.begin();i != swaps.end();i++) {
		//printf("SWAP BETWEEN %d and %d\n", i->first, i->second);
	}
	*/

	/*Test for lazy circuit builder*/
	betterLazySwapCircuitBuilder(architectureEdges, firstGates, distanceMatrix);
	//Exit the program:
	return 0;
}

// another implementation for a mapper for the initial mapping project
	//just returns a single perfect mapping if found
map<int, int> perfectMapper(map<int, set<int>> architectureEdges, set<GateNode*> startSet) {
	//we first crop the edge graph and proceed as in the partial mapper 
	//however we will use dfs for this mapper
	map<int,set<int>> copyArchitectureEdges =architectureEdges;
	
	map<int,set<int>> circuitEdges;
	
	set<GateNode*> closedSet;
	//this is a vector to insert accounted gates into, will come into play when we actually start mapping
	vector<GateNode*> mappingQueue;

		
	queue<GateNode*> startGates;
	set<GateNode*> start = startSet;
	while (start.size() > 0) {
		////printf("current size: %d\n",start.size());
		for (auto i = start.begin();i!=start.end();) {
			
			if (start.find((*i)->controlParent) == start.end() && start.find((*i)->targetParent) == start.end()) {
				//then we add it into queue to get a level order
				startGates.push((*i));
				i=start.erase(i);
				////printf("gate pushed!\n");	
			} else{
				i++;
			}
		}
		////printf("exited loop!\n");
	}

	//going gate by gate through each level, adding edges and then checking if we went too far
	//big note here: gates in each level do not depend on each other, so we know we can always map level by level
	while (startGates.size() > 0) {
		//adding the current level of the dependence graph to our circuit graph and updating node degrees

		GateNode* chosen = startGates.front();
		startGates.pop();
		////printf("Chosen gate: %d control and %d target\n",chosen->control,chosen->target);
		if (closedSet.find(chosen) != closedSet.end()) {
			//then we already saw this
			continue;
		}
		//we will not consider 1 qubit gates here
		if (chosen->control==-1){
			//still might need to enqueue target child
				//we consider single qubit gates as 
					//"automatically satisfied"
			closedSet.insert(chosen);
			mappingQueue.push_back(chosen);
			if (chosen->controlChild!=NULL){
				////printf("ITS A TRAP!\n");
			}
			if (chosen->targetChild!=NULL){
				GateNode* child = chosen->targetChild;
				if (child->control == chosen->target){
					if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end()){
						startGates.push(child);
					}
				} else{
					if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end()){
						startGates.push(child);
					}
				}
				
			}
			
			
			continue;
		}
		
		////printf("starting add edge\n");
		circuitEdges = addGateEdge(chosen, circuitEdges);
		unsigned int circuitMinDegree=circuitEdges.size();
		unsigned int circuitMaxDegree=0;
		for (auto i=circuitEdges.begin();i!=circuitEdges.end();i++){
			if ((i->second).size()>circuitMaxDegree){
				circuitMaxDegree = (i->second).size();
			}
			if ((i->second).size()<circuitMinDegree){
				circuitMinDegree = (i->second).size();
			}
		}
		////printf("starting prune edges\n");
		pair<queue<pair<int,int>>,map<int,set<int>>> pruneRes = pruneArchitectureEdges(copyArchitectureEdges,circuitMinDegree);
		queue<pair<int,int>> changes = pruneRes.first;
		copyArchitectureEdges = pruneRes.second;

		////printf("passed gateEdgeAdd and prune!\n");
		////printf("sanity check: \n");	
		////printf("min degree of circ is %d\n",circuitMinDegree);
		unsigned int archMinDegree=copyArchitectureEdges.size();
		for (auto i =copyArchitectureEdges.begin() ; i!=copyArchitectureEdges.end();i++){
			if ((i->second).size() < archMinDegree){
				archMinDegree = (i->second).size();
			}
		}
		////printf("min degree of arch is %d\n",archMinDegree);
		
		
		//suppose we went too far
			//then we rollback and continue with the queue
			// keep doing this until queue is emtpy --> that is our maximal partition
		
		//getting max degree of the architecture nodes
		unsigned int architectureMaxDegree = 0;
		for (auto i = copyArchitectureEdges.begin();i != copyArchitectureEdges.end();i++) {
			if ((i->second).size() > architectureMaxDegree) {
				architectureMaxDegree = (i->second).size();
			}
		}
		if (wentTooFar(circuitEdges,copyArchitectureEdges,circuitMaxDegree,architectureMaxDegree)) {
			////printf("GONE TOO FAR! MEANS NO PERFECT INITIAL MAPPING!\n");
			//returning empty mapping
			return map<int,int>();
		}
		else {
			//then we can update (add children of the block to the queue
			//adding the chosen gate to the mapping queue
			mappingQueue.push_back(chosen);
			closedSet.insert(chosen);
			if (chosen->controlChild!=NULL){
				GateNode* child = chosen->controlChild;
				if (child->control == chosen->control){
					//then we need to check target parent
					if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end()){
						startGates.push(child);
					}
				} else {
					//then we need to check control parent
					if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end()){
						startGates.push(child);
					}
				}
			}

			if (chosen->targetChild!=NULL){
				GateNode* child = chosen->targetChild;
				
				if (child->control == chosen->target){
					//then we need to check target parent
					if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end()){
						startGates.push(child);
					}
				} else {
					
					//then we need to check control parent
					if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end()){
						startGates.push(child);
					}
				}
			}
			
			
		}
	}

	//now mapping queue is a vector of gates to map in order (if we reached here, we might be able to map every vertex)
	stack<tuple<map<int, int>,set<int>, int>> mappingStack;
	//stack of mappings and counter to which gatenode in the vector was satisfied, once we reach a mapping with all gate nodes satisfied, we can just return it
	// if the stack becomes empty, then we can return NULL

	//pushing initial setup onto stack (for gate 0) 
	int count = 0;
	// we do not want to map for single qubit gates
	while (mappingQueue[count]->control == -1) {
		count++;
	}
	//now count points to a gate in mappingQueue that is a 2 qubit gate 
	//we push possible pairings of qubits to satisfy it into the stack
	int target = mappingQueue[count]->target;
	int control = mappingQueue[count]->control;
	
	for (auto i = architectureEdges.begin();i != architectureEdges.end();i++) {
		if ((i->second).size() >= (circuitEdges[control]).size()) {
				//this is a possible control
			for (auto j = (i->second).begin();j != (i->second).end();j++) {
				if (architectureEdges[*j].size() >= (circuitEdges[target]).size()) {
					//this is a possible target
					// we will push the entry onto the stack
					map<int, int> currMapping;
					currMapping[control] = i->first;
					currMapping[target] = *j;

					set<int> mappedArch;
					mappedArch.insert(i->first);
					mappedArch.insert(*j);

					int gateNum = count;
					mappingStack.push(make_tuple(currMapping, mappedArch,gateNum));
				}
			}
		}
	}

	while (mappingStack.size() > 0) {
		//getting the top mapping
		tuple<map<int, int>,set<int>, int> stackEntry = mappingStack.top();
		mappingStack.pop();
		map<int, int> currMapping = get<0>(stackEntry);
		set<int> mappedArchQubits = get<1>(stackEntry);
		unsigned int gateReached = get<2>(stackEntry);
		if (gateReached == mappingQueue.size() - 1) {
			//then we got a perfect initial mapping
			return currMapping;
		}

		GateNode* toSatisfy = mappingQueue[gateReached + 1];

		if (toSatisfy->control == -1) {
			//single qubit gate
			//we can just increment counter and skip
			mappingStack.push(make_tuple(currMapping, mappedArchQubits, gateReached + 1));
			continue;
		}

		// we have 3 cases now
			//1) both qubits are already mapped 
			//2) one of the qubits are mapped
			//3) none of the qubits are mapped
		if (currMapping.find(toSatisfy->control) != currMapping.end() && currMapping.find(toSatisfy->target) != currMapping.end()) {
				////printf("case 1\n");
				//then both are already mapped, we just check if they satisfy the new gate to satisfy
				//namely the architecture nodes they are mapped to should have an edge between them
				set<int> edges = copyArchitectureEdges[currMapping[toSatisfy->control]];
				if (edges.find(currMapping[toSatisfy->target])!=edges.end()){
					//then the edge is there
					mappingStack.push(make_tuple(currMapping,mappedArchQubits, gateReached+1));
				} 				
			}
			else if (currMapping.find(toSatisfy->control) != currMapping.end() || currMapping.find(toSatisfy->target) != currMapping.end()) {
				////printf("case 2\n");
				//then only 1 of the two qubits in the gate are mapped
				//we find possible second qubits for the mapping to complete the gate
					//if none found, then mapping isnt propagated
				int mappedQubit = 0;
				int unmappedQubit = 0;
				if (currMapping.find(toSatisfy->control) != currMapping.end()) {
					mappedQubit = toSatisfy->control;
					unmappedQubit = toSatisfy->target;
				}
				else {
					mappedQubit = toSatisfy->target;
					unmappedQubit = toSatisfy->control;
				}

				//finding possible second mappings which will satisfy the gate
				set<int> edges = copyArchitectureEdges[currMapping[mappedQubit]];
				for (auto i= edges.begin();i!=edges.end();i++){
					if (copyArchitectureEdges[*i].size()>= circuitEdges[unmappedQubit].size() && mappedArchQubits.find(*i)==mappedArchQubits.end()){
						//mapping is possible here
						map<int, int> newMapping = currMapping;
						set<int> archMapping = mappedArchQubits;
						archMapping.insert(*i);
						newMapping.insert(make_pair(unmappedQubit, *i));
						mappingStack.push(make_tuple(newMapping,archMapping, gateReached+1));
					}
				}
			}
			else {
				//then none of the 2 qubits in the gate are mapped
				////printf("case 3\n");
					//we need to extend this mapping by finding possible pairs of qubits
					//so we will first find eligible controls, and for each control an eligible target
				for (auto i = copyArchitectureEdges.begin();i!=copyArchitectureEdges.end();i++){
					if (mappedArchQubits.find(i->first)==mappedArchQubits.end() && (i->second).size()>=(circuitEdges[toSatisfy->control]).size()){
					//possible pairing for the control, we need to also go through possible pairings for the target
						set<int> edges = i->second;
						for (auto j=edges.begin();j!=edges.end();j++){
							if (mappedArchQubits.find(*j)==mappedArchQubits.end() && (copyArchitectureEdges[*j]).size()>=(circuitEdges[toSatisfy->target]).size()  ){
								//this is a good pairing, we should propagate
								map<int, int> newMapping = currMapping;
								set<int> newArchMapping = mappedArchQubits;
								newMapping.insert(make_pair(toSatisfy->control, i->first));
								newMapping.insert(make_pair(toSatisfy->target, *j));
								newArchMapping.insert(i->first);
								newArchMapping.insert(*j);
								mappingStack.push(make_tuple(newMapping, newArchMapping, gateReached+1));

							}
						}
					}
				}
				
				
			}

		}	

		//if we reached here, then no perfect mapping exists
		return map<int,int>();

}

//maximal mapper that returns matchings for the next maximal portion of the dependence graph
	//returns a pair
		//the first element is a list of all maximal matchings
		//the second element is jsut a set of mapped architecture nodes --> for quick lookup
		//the third element is the queue of gatenodes which could not be matched and should be considered the "start" for the next iteration
		//this is returned so we can easily start the next iteration of mapping
vector <tuple<map<int, int>, set<int>, set<GateNode*>> > maximalMapper( map<int, set<int>> architectureEdges, set<GateNode*> startSet) {
	//we first grab a copy of the coupling graph and associated nodeDegree vector
	////printf("entered maximal mapper!\n");
	//multimap<int, int> copyArchitectureEdges = architectureEdges;
	//smart pointer copy of architectureEdges
	map<int,set<int>> copyArchitectureEdges =architectureEdges;
	//map<int,int> copyArchitectureNodeDegrees = architectureNodeDegrees;

	//allocating graph structure to add level by level in the dependence graph
	//multimap<int, int> circuitEdges ;
	//smart pointer to circuit graph
	map<int,set<int>> circuitEdges;
	//map<int,int> circuitNodeDegrees ;


	//keeping track of considered gates (we do not want to consider gates twice)
	set<GateNode*> closedSet;
	//this is a queue to insert accounted gates into, will come into play when we actually start mapping
	deque<GateNode*> mappingQueue;
	//vector<GateNode*> mappingQueue;
	//set<GateNode*> mappingSet;
	set<GateNode*> couldntMap;

		
	queue<GateNode*> startGates;
	set<GateNode*> start = startSet;
	////printf("input gates size: %d",startSet.size());
	////printf("pushing start gates in order into queue!\n");
	while (start.size() > 0) {
		////printf("current size: %d\n",start.size());
		for (auto i = start.begin();i!=start.end();) {
			
			if (start.find((*i)->controlParent) == start.end() && start.find((*i)->targetParent) == start.end()) {
				//then we add it into queue to get a level order
				startGates.push((*i));
				i=start.erase(i);
				////printf("gate pushed!\n");	
			} else{
				i++;
			}
		}
		////printf("exited loop!\n");
	}

	//going gate by gate through each level, adding edges and then checking if we went too far
	//big note here: gates in each level do not depend on each other, so we know we can always map level by level
	////printf("Determining maximal possible mapping of gates!\n");
	////printf("Starting size of startGates: %d!\n",startGates.size());
	while (startGates.size() > 0) {
		//adding the current level of the dependence graph to our circuit graph and updating node degrees

		GateNode* chosen = startGates.front();
		startGates.pop();
		////printf("Chosen gate: %d control and %d target\n",chosen->control,chosen->target);
		if (chosen->controlChild!=NULL && chosen->controlChild->control==8 && chosen->controlChild->target==12){
			////printf("THIS GATE IS THE CULPRIT!\n");
		}
		if (chosen->targetChild!=NULL && chosen->targetChild->control==8 && chosen->targetChild->target==12){
			////printf("THIS GATE IS THE CULPRIT!\n");
		}
		if (closedSet.find(chosen) != closedSet.end()) {
			//then we already saw this
			continue;
		}
		//we will not consider 1 qubit gates here
		if (chosen->control==-1){
			//still might need to enqueue target child
				//we consider single qubit gates as 
					//"automatically satisfied"
			closedSet.insert(chosen);
			mappingQueue.push_back(chosen);
			if (chosen->controlChild!=NULL){
				////printf("ITS A TRAP!\n");
			}
			if (chosen->targetChild!=NULL){
				GateNode* child = chosen->targetChild;
				if (child->control == chosen->target){
					if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end() || ( startSet.find(child->targetParent)==startSet.end())){
						startGates.push(child);
					}
				} else{
					if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end() || (startSet.find(child->controlParent)==startSet.end())){
						startGates.push(child);
					}
				}
				
			}
			
			
			continue;
		}

		//going to closed set and seeing if there are any target, control repeats
			//if so, they are automatically satisfied
		bool repeat = false;
		for (auto z = closedSet.begin();z != closedSet.end();z++) {
			if ((chosen->control == (*z)->control && chosen->target == (*z)->target) || (chosen->target == (*z)->control && chosen->control == (*z)->target)) {
				//this gate is a "repeat" gate in terms of graph cutting, we can just consider it satisfied and add its children
				mappingQueue.push_back(chosen);
				closedSet.insert(chosen);
				if (chosen->controlChild!=NULL){
					GateNode* child = chosen->controlChild;
					if (child->control == chosen->control){
						//then we need to check target parent
						if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end() || startSet.find(child->targetParent)==startSet.end()){
							startGates.push(child);
						}
					} else {
						//then we need to check control parent
						if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end() || startSet.find(child->controlParent)==startSet.end()){
							startGates.push(child);
						}
					}
				}

				if (chosen->targetChild!=NULL && chosen->targetChild != chosen->controlChild){
					GateNode* child = chosen->targetChild;
					
					if (child->control == chosen->target){
						//then we need to check target parent
						if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end() || startSet.find(child->targetParent)==startSet.end()){
							startGates.push(child);
						}
					} else {
						
						//then we need to check control parent
						if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end() || startSet.find(child->controlParent)==startSet.end()){
							startGates.push(child);
						}
					}
				}
				repeat = true;
				break;
			}
			
		}
		if (repeat) {
			continue;
		}
		
		////printf("starting add edge\n");
		circuitEdges = addGateEdge(chosen, circuitEdges);
		unsigned int circuitMinDegree=circuitEdges.size();
		unsigned int circuitMaxDegree=0;
		for (auto i=circuitEdges.begin();i!=circuitEdges.end();i++){
			if ((i->second).size()>circuitMaxDegree){
				circuitMaxDegree = (i->second).size();
			}
			if ((i->second).size()<circuitMinDegree){
				circuitMinDegree = (i->second).size();
			}
		}
		////printf("starting prune edges\n");
		pair<queue<pair<int,int>>,map<int,set<int>>> pruneRes = pruneArchitectureEdges(copyArchitectureEdges,circuitMinDegree);
		queue<pair<int,int>> changes = pruneRes.first;
		copyArchitectureEdges = pruneRes.second;

		////printf("passed gateEdgeAdd and prune!\n");
		////printf("sanity check: \n");	
		////printf("min degree of circ is %d\n",circuitMinDegree);
		unsigned int archMinDegree=copyArchitectureEdges.size();
		for (auto i =copyArchitectureEdges.begin() ; i!=copyArchitectureEdges.end();i++){
			if ((i->second).size() < archMinDegree){
				archMinDegree = (i->second).size();
			}
		}
		////printf("min degree of arch is %d\n",archMinDegree);
		
		
		//suppose we went too far
			//then we rollback and continue with the queue
			// keep doing this until queue is emtpy --> that is our maximal partition
		
		//getting max degree of the architecture nodes
		unsigned int architectureMaxDegree = 0;
		for (auto i = copyArchitectureEdges.begin();i != copyArchitectureEdges.end();i++) {
			if ((i->second).size() > architectureMaxDegree) {
				architectureMaxDegree = (i->second).size();
			}
		}
		if (wentTooFar(circuitEdges,copyArchitectureEdges,circuitMaxDegree,architectureMaxDegree)) {
			////printf("GONE TOO FAR!\n");
			//rollback 
			copyArchitectureEdges=rollBack(changes, copyArchitectureEdges);
			//undoing circuit node add and edges/degree
			circuitEdges[chosen->control].erase(chosen->target);
			circuitEdges[chosen->target].erase(chosen->control);
			if (circuitEdges[chosen->control].size() == 0) {
				circuitEdges.erase(chosen->control);
			}

			if (circuitEdges[chosen->target].size() == 0) {
				circuitEdges.erase(chosen->target);
			}

			//we couldnt map the given gate, so we will push it to a list of nodes 
			couldntMap.insert(chosen);
			//we do not want to consider this gate again
			closedSet.insert(chosen);
		}
		else {
			//then we can update (add children of the block to the queue
			//adding the chosen gate to the mapping queue
			mappingQueue.push_back(chosen);
			closedSet.insert(chosen);
			if (chosen->controlChild!=NULL){
				GateNode* child = chosen->controlChild;
				if (child->control == chosen->control){
					//then we need to check target parent
					if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end() || startSet.find(child->targetParent)==startSet.end()){
						startGates.push(child);
					}
				} else {
					//then we need to check control parent
					if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end() || startSet.find(child->controlParent)==startSet.end()){
						startGates.push(child);
					}
				}
			}

			if (chosen->targetChild!=NULL && chosen->targetChild != chosen->controlChild){
				GateNode* child = chosen->targetChild;
				
				if (child->control == chosen->target){
					//then we need to check target parent
					if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end() || startSet.find(child->targetParent)==startSet.end()){
						startGates.push(child);
					}
				} else {
					
					//then we need to check control parent
					if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end() || startSet.find(child->controlParent)==startSet.end()){
						startGates.push(child);
					}
				}
			}
			
			
		}
	}
	

	//after the queue is empty, we can start mapping and now we will get a maximal mapping for this level
		//we map gate by gate in the order of enqueuing in the mappingQueue object
	//stack is of the form ((logical q_1, physical mapping) , (logical q_2, physical mapping))
		//this is because we are setting the maximal mapping gate by gate
		//therefore, we need to consider pairings of both the target and the control 
			//if target and control is unmapped
			//if only target is unmapped
			// if only control is unmapped
			//if both control and target are mapped
	//debug: maximalMappingQueue order check	
	int count=0;
	int max = mappingQueue.size();
	set<GateNode*> debugSet;
	while (count!=max){
		GateNode* forMapping = mappingQueue.front();
		mappingQueue.pop_front();
		mappingQueue.push_back(forMapping);
		if ((forMapping->controlParent !=NULL && debugSet.find(forMapping->controlParent)==debugSet.end()) || (forMapping->targetParent!=NULL && debugSet.find(forMapping->targetParent)==debugSet.end())){
			////printf("ORDER NOT SATISFIED! for control %d and target %d\n",forMapping->control,forMapping->target);
			if ((forMapping->controlParent !=NULL && forMapping->controlParent->control==-1) ||(forMapping->targetParent!=NULL && forMapping->targetParent->control==-1) ){
				////printf("AS I SUSPECTED!\n");
			}
		}
		debugSet.insert(forMapping);
		count++;
	}

	////printf("CONFIRMING THE PRUNE!\n");
	////printf("HARDWARE EDGES!!!!\n");
	for (auto z = copyArchitectureEdges.begin();z != copyArchitectureEdges.end();z++) {
		for (auto p = z->second.begin();p != z->second.end();p++) {
			////printf("WE HAVE EDGE FROM HARDWARE QUBIT %d to HARDWARE QUBIT %d\n", z->first, *p);
		}
	}

	////printf("CIRCUIT EDGES!!!!\n");
	for (auto z = circuitEdges.begin();z != circuitEdges.end();z++) {
		for (auto p = z->second.begin();p != z->second.end();p++) {
			//printf("WE HAVE EDGE FROM LOGICAL QUBIT %d to LOGICAL QUBIT %d\n", z->first, *p);
		}
	}


	//our stack is of the form <(control, pairing),(target,pairing)>
	//stack < pair<pair<int, int>,pair<int,int>>> maximalMappingStack;
	//using a queue to go gate by gate and extend our maximal mappings
	//printf("starting actual isomorphism finding!\n");
	queue < tuple<map<int, int>,set<int>, set<GateNode*>>> maximalMappingQueue;
	//keeping track of where we are in the mapping
	set<GateNode*> otherClosedSet; 
	
	//we do the first iteration
	GateNode* firstGate = mappingQueue.front();	
	mappingQueue.pop_front();
	otherClosedSet.insert(firstGate);
	int control = firstGate->control;
	int target = firstGate->target;
	//printf("now we hit gate with control %d and target %d\n",control,target);
	if (control==-1){
		//then we will just push an empty mapping onto the queue
		map<int,int> emptyMapping;
		set<int> emptyArch;
		set<GateNode*> initial;
		initial.insert(firstGate);
		maximalMappingQueue.push(make_tuple(emptyMapping,emptyArch, initial));

	} else{
		auto iter = copyArchitectureEdges.begin();
		while (iter != copyArchitectureEdges.end()) {
			if ((iter->second).size() >= circuitEdges[control].size()) {
				//then we found a physical qubit that can accomodate the logical qubit "control" in this maximal mapping
				//we need to check the neighbors of this qubit for possible "target" mappings
				for (auto i = (iter->second).begin();i!=(iter->second).end();i++){
					if (copyArchitectureEdges[*i].size()>= circuitEdges[target].size()){
						//printf("REACHED VERY INNER PART!\n");
						//then this is a possible mapping for the target logical qubit
						map<int, int> suggestedInitialMapping;
						set<int> mappedArchitectureQubits;
						suggestedInitialMapping.insert(make_pair(control, iter->first));
						mappedArchitectureQubits.insert(iter->first);
						suggestedInitialMapping.insert(make_pair(target, *i));
						mappedArchitectureQubits.insert(*i);
						set<GateNode*> currentlySatisfied;
						currentlySatisfied.insert(firstGate);
						//pushing the mapping and satisfied gate to the queue 	
						maximalMappingQueue.push(make_tuple(suggestedInitialMapping,mappedArchitectureQubits, currentlySatisfied));

					}
				}
				
			}
			iter++;
		}
	}
	//printf("INITIAL QUEUE SIZE: %lu\n", maximalMappingQueue.size());
	//starting mapping loop
	//keeping track of the best mappings at end of every step
	unsigned int lastMapped =1;
	while (mappingQueue.size() > 0) {
		GateNode* toSatisfy = mappingQueue.front();
		mappingQueue.pop_front();
		//printf("toSatisfy: %d control and %d target\n",toSatisfy->control,toSatisfy->target);
		 	

		//we go through all mappings in the queue and see if we can propogate them to the next gate
			//we can easily check if parent gates are satisfied in the current mapping or not with the extra set given
		//queue < tuple<map<int, int>, set<int>,set<GateNode*>>> resultMaximalMappingQueue;
		unsigned int internalStep =0;
		unsigned int size = maximalMappingQueue.size();
		//debug
		//printf("queue size: %d\n",size);
		while (size!=0) {
			tuple<map<int, int>, set<int>, set<GateNode*>> queueEntry = maximalMappingQueue.front();
			maximalMappingQueue.pop();
			size--;
			
			//seeing if we can propagate the mapping to the next gate, and inserting it in the queue
			map<int, int> currMapping = get<0>(queueEntry);
			set<int> mappedArchQubits = get<1>(queueEntry);
			set<GateNode*> satisfiedGates = get<2>(queueEntry);
			//need to throw out mapping if some conditions are not satisfied
			//perfect mapping condition for now
			if (satisfiedGates.size() <lastMapped){
				//then we can safely throw this out for now because there are better mappings
				
				continue;
			}

			
			
			
			if ((debugSet.find(toSatisfy->controlParent)!= debugSet.end() && satisfiedGates.find(toSatisfy->controlParent)==satisfiedGates.end()) || (debugSet.find(toSatisfy->targetParent)!=debugSet.end() && satisfiedGates.find(toSatisfy->targetParent)==satisfiedGates.end())){
				
				//parents not satisfied
				maximalMappingQueue.push(queueEntry);
				if ((debugSet.find(toSatisfy->controlParent)!= debugSet.end() && satisfiedGates.find(toSatisfy->controlParent)==satisfiedGates.end())){
					//printf("unsatisfied parent: %d control and %d target\n",toSatisfy->controlParent->control,toSatisfy->controlParent->target);
				}
				if ((debugSet.find(toSatisfy->targetParent)!= debugSet.end() && satisfiedGates.find(toSatisfy->targetParent)==satisfiedGates.end())){
					//printf("unsatisfied parent: %d control and %d target\n",toSatisfy->targetParent->control,toSatisfy->targetParent->target);
				}
				continue;
			}
			//checking if we got a single qubit gate
				//if so, we automatically consider it satisfied
			if (toSatisfy->control==-1){
				//automatically satisfied single qubit gate
				////printf("AUTOMATICALLY SATISFIED!\n");
				set<GateNode*> newSatisfied = satisfiedGates;	
				newSatisfied.insert(toSatisfy);
				if (newSatisfied.size()>internalStep){
					internalStep = newSatisfied.size();
				}
				maximalMappingQueue.push(make_tuple(currMapping,mappedArchQubits,newSatisfied));
				continue;
			}
			

			//dependencies are resolved if we reached here
				//we need to check 3 cases
				//1) both qubits are mapped
				//2) only 1 of the two qubits are mapped
				//3) none of the qubits are mapped
			if (currMapping.find(toSatisfy->control) != currMapping.end() && currMapping.find(toSatisfy->target) != currMapping.end()) {
				////printf("case 1\n");
				//then both are already mapped, we just check if they satisfy the new gate to satisfy
				//namely the architecture nodes they are mapped to should have an edge between them
				set<int> edges = copyArchitectureEdges[currMapping[toSatisfy->control]];
				if (edges.find(currMapping[toSatisfy->target])!=edges.end()){
					//then the edge is there
					satisfiedGates.insert(toSatisfy);
					if (satisfiedGates.size()>internalStep){
						internalStep = satisfiedGates.size();
					}
					maximalMappingQueue.push(make_tuple(currMapping,mappedArchQubits, satisfiedGates));
				} else{
					//the edge is not present
					if (satisfiedGates.size()>internalStep){
						internalStep = satisfiedGates.size();
					}

					maximalMappingQueue.push(queueEntry);
				}
				
			}
			else if (currMapping.find(toSatisfy->control) != currMapping.end() || currMapping.find(toSatisfy->target) != currMapping.end()) {
				////printf("case 2\n");
				//then only 1 of the two qubits in the gate are mapped
				//we find possible second qubits for the mapping to complete the gate
					//if none found, then mapping isnt propagated
				int mappedQubit = 0;
				int unmappedQubit = 0;
				if (currMapping.find(toSatisfy->control) != currMapping.end()) {
					mappedQubit = toSatisfy->control;
					unmappedQubit = toSatisfy->target;
				}
				else {
					mappedQubit = toSatisfy->target;
					unmappedQubit = toSatisfy->control;
				}

				//finding possible second mappings which will satisfy the gate
				bool didPropagate = false;
				set<int> edges = copyArchitectureEdges[currMapping[mappedQubit]];
				for (auto i= edges.begin();i!=edges.end();i++){
					if (copyArchitectureEdges[*i].size()>= circuitEdges[unmappedQubit].size() && mappedArchQubits.find(*i)==mappedArchQubits.end()){
						//mapping is possible here
						didPropagate = true;
						map<int, int> newMapping = currMapping;
						set<int> archMapping = mappedArchQubits;
						archMapping.insert(*i);
						newMapping.insert(make_pair(unmappedQubit, *i));
						set<GateNode*> newSatisfied = satisfiedGates;
						newSatisfied.insert(toSatisfy);
						if (newSatisfied.size()>internalStep){
							internalStep = newSatisfied.size();
						}
						maximalMappingQueue.push(make_tuple(newMapping,archMapping, newSatisfied));
					}
				}
				
				if (!didPropagate) {
					//moving the old entry in case we couldnt propogate it to the next gate
					if (satisfiedGates.size()>internalStep){
						internalStep = satisfiedGates.size();
					}
					maximalMappingQueue.push(queueEntry);
				}
			}
			else {
				//then none of the 2 qubits in the gate are mapped
				////printf("case 3\n");
					//we need to extend this mapping by finding possible pairs of qubits
					//so we will first find eligible controls, and for each control an eligible target
				bool propagated=false;
				for (auto i = copyArchitectureEdges.begin();i!=copyArchitectureEdges.end();i++){
					if (mappedArchQubits.find(i->first)==mappedArchQubits.end() && (i->second).size()>=(circuitEdges[toSatisfy->control]).size()){
					//possible pairing for the control, we need to also go through possible pairings for the target
						set<int> edges = i->second;
						for (auto j=edges.begin();j!=edges.end();j++){
							if (mappedArchQubits.find(*j)==mappedArchQubits.end() && (copyArchitectureEdges[*j]).size()>=(circuitEdges[toSatisfy->target]).size()  ){
								//this is a good pairing, we should propagate
								propagated = true;
								map<int, int> newMapping = currMapping;
								set<int> newArchMapping = mappedArchQubits;
								set<GateNode*> newSatisfied = satisfiedGates;
								newMapping.insert(make_pair(toSatisfy->control, i->first));
								newMapping.insert(make_pair(toSatisfy->target, *j));
								newArchMapping.insert(i->first);
								newArchMapping.insert(*j);
								newSatisfied.insert(toSatisfy);
									if (newSatisfied.size()>internalStep){
					internalStep = newSatisfied.size();
    }
								maximalMappingQueue.push(make_tuple(newMapping, newArchMapping, newSatisfied));

							}
						}
					}
				}
				
				if (!propagated) {
					//pushing old status then
					if (satisfiedGates.size()>internalStep){
						internalStep = satisfiedGates.size();
					}
					maximalMappingQueue.push(queueEntry);
				}
			}

		}
		
		//setting our new queue to be the result queue (we basically just pushed propagations to a new queue and now we are setting it as the old queue)
		//maximalMappingQueue = resultMaximalMappingQueue;
		otherClosedSet.insert(toSatisfy);
		//updating counter of number of mappings in best mapping so far
		lastMapped = internalStep;
	}
	//printf("PRINTING GATES CONSIDERED IN THIS ROUND OF MAPPING:\n");
	for (auto z = otherClosedSet.begin();z != otherClosedSet.end();z++) {
		//cout << "Gate name: " + (*z)->name << " control: " << (*z)->control << " target: " << (*z)->target <<"\n";
	}
	//printf("FINISHED MAPPING DOING PREP NOW!\n");
	vector <tuple<map<int, int>, set<int>, set<GateNode*>> > results;
	queue <tuple<map<int, int>, set<int>, set<GateNode*>>> resultQueue;
	unsigned int mostMatched = 0;
	while (maximalMappingQueue.size() > 0) {
		tuple<map<int, int>, set<int>, set<GateNode*>> queueEntry = maximalMappingQueue.front();
		maximalMappingQueue.pop();
		if (get<2>(queueEntry).size() > mostMatched) {
			mostMatched = get<2>(queueEntry).size();
		}
		resultQueue.push(queueEntry);
		
		
	}
	int numGatesMapped = 0;
	while (resultQueue.size() > 0) {
		tuple<map<int, int>, set<int>, set<GateNode*>> queueEntry = resultQueue.front();
		resultQueue.pop();
		if (get<2>(queueEntry).size() < mostMatched) {
			continue;
		}
		
		set<GateNode*> completedGates = get<2>(queueEntry);
		if (completedGates.size() > numGatesMapped) {
			numGatesMapped = completedGates.size();
		}
		set<GateNode*> finalSet;
		for (auto j = otherClosedSet.begin();j != otherClosedSet.end();j++) {
			if (completedGates.find((*j)) == completedGates.end()) {
				//then this wasnt mapped
				finalSet.insert(*j);
			}
		}

		for (auto j = couldntMap.begin();j != couldntMap.end();j++) {
			finalSet.insert(*j);
		}
		results.push_back(make_tuple(get<0>(queueEntry), get<1>(queueEntry), finalSet));

	}
		
	//now our maximalMappingQueue holds all the generated partial mappings, we have to do some cleanup to get the best ones
		//in addition, we need to do some setup for the future initial gates to consider (because we will run this method multiple times)
	//returning the best partial mappings, and the gates for the next iteration are already modified in firstGates
	//printf("RETURNING RESULTS!\n");
	cout << "Gates mapped: " << numGatesMapped << "\n";
	cout << "COULDNT MAP SIZE! : " << couldntMap.size() << " \n";
	cout << "Gates Considered Size: " << otherClosedSet.size() << "\n";
	return results;
	
}

//helper function to get the next level of the dependence graph, given the previous level
set<GateNode*> getNextLevel(set<GateNode*> previousLevel) {
	// by definition, gates in a level of a dependence graph do not depend on each other
	// therefore, we get all the children of the previousLevel
		//if any children are children of other gates in the set, we can remove them
	set<GateNode*> toReturn ;
	auto iter = previousLevel.begin();
	while (iter != previousLevel.end()) {
		GateNode* prevGate = *iter;
		//adding children to the return set
		if (prevGate->controlChild != NULL) {
			toReturn.insert(prevGate->controlChild);
		}

		if (prevGate->targetChild != NULL) {
			toReturn.insert(prevGate->targetChild);
		}
		iter++;
	}

	//removing gates which have parents in the level set
	set<GateNode*> toRemove;
	iter = previousLevel.begin();
	while (iter!=previousLevel.end()){
		GateNode* toCheck = *iter;
		if (toReturn.find(toCheck->controlParent) != toReturn.end() || toReturn.find(toCheck->targetParent) != toReturn.end()) {
			//adding this gate to list of gates to remove
			toRemove.insert(toCheck);
		}
		iter++;
	}

	//removing all the gates that have parents in the set already (those need to be satisfied first)
	iter = toRemove.begin();
	while (iter != toRemove.end()) {
		toReturn.erase(*iter);
		iter++;
	}
	//returning set of next level of gates in dependence graph
	return toReturn;
}


//helper method to figure out if we went too far in our initial mapping or not
	//returns true if we need to roll back
	//false otherwise (as per name "did I go to far or not?")
bool wentTooFar(map<int,set<int>> circuitEdges, map<int,set<int>>architectureEdges,int circuitMaxDegree, int architectureMaxDegree) {
	//another prior check for going through both lists with two "fingers" and removing elements as I go for mappings
	//if I cannot remove all elements from the copyOfCircuitDegrees vector, then no complete initial mapping exists
	//we will need to roll back changes 
	//now we have left architecture nodes with at least the degree needed to satisfy the circuit requirements
	if (architectureEdges.size() < circuitEdges.size()) {
		//then we do not have the required number of nodes to try and find a perfect initial mapping
		cout << "NO MAPPING POSSIBLE! SPOT 1\n";
		return true;
	}

	if (circuitMaxDegree > architectureMaxDegree) {
		//then we cut out too much from the architecture graph
		cout << "NO MAPPING POSSIBLE ! SPOT 2\n";
		return true;
	}

	//sorting copies of both lists from least to greatest node degrees
		//we only care about the values, so we will use vector<int>
	vector<int> copyOfArchDegrees;
	vector<int> copyOfCircuitDegrees;
	for (auto i = architectureEdges.begin();i != architectureEdges.end();i++) {
		copyOfArchDegrees.push_back((i->second).size());
	}
	for (auto i = circuitEdges.begin();i != circuitEdges.end();i++) {
		copyOfCircuitDegrees.push_back((i->second).size());
	}
	sort(copyOfArchDegrees.begin(), copyOfArchDegrees.end());
	sort(copyOfCircuitDegrees.begin(), copyOfCircuitDegrees.end());
	
	auto architectureNodesFinger = copyOfArchDegrees.begin();
	auto circuitNodesFinger = copyOfCircuitDegrees.begin();

	while (architectureNodesFinger != copyOfArchDegrees.end() &&
		circuitNodesFinger != copyOfCircuitDegrees.end()) {
		if ((*circuitNodesFinger) <= (*architectureNodesFinger)) {
			//advancing both iterators
			circuitNodesFinger++;
			architectureNodesFinger++;
		}
		else {
			//only advancing the architectureNodesFinger
			architectureNodesFinger++;
		}
	}

	//there is not possible complete initial mapping if the circuitNodesFinger did not reach the end
	if (circuitNodesFinger != copyOfCircuitDegrees.end()) {
		//cout << "NO POSSIBLE MATCHING! SPOT 3\n";
		return true;
	}
	return false;
}

//adds an edge to our maximal mapping graph and returns the smallest and largest degree accordingly
map<int,set<int>> addGateEdge(GateNode* chosen, map<int,set<int>>circuitEdges) {
	int target = chosen->target;
	int control = chosen->control;
	if (circuitEdges.find(target)==circuitEdges.end()){
		set<int> edges;
		edges.insert(control);
		circuitEdges[target]=edges;
	} else{
		(circuitEdges[target]).insert(control);
	}

	if (circuitEdges.find(control)==circuitEdges.end()){
		set<int> edges;
		edges.insert(target);
		circuitEdges[control]=edges;
	} else{
		(circuitEdges[control]).insert(target);
	}
	

	
	return circuitEdges;
}

//given a max and min degree of the circuit graph I have defined, we can prune the coupling graph
	//returns all the changes made to the graph (in case we need to roll back)
 pair<queue<pair<int,int>>,map<int,set<int>>> pruneArchitectureEdges(map<int, set<int>>architectureEdges, unsigned int circuitMinDegree) {
	//iterating through architectureNodeDegrees and removing any vertex from the graph that is below the min degree	
	 //recording list of changes made
	 queue<pair<int, int>> changes;
	bool changeMade = false;
	do{

		auto iter = architectureEdges.begin();
		changeMade=false;
		while (iter != architectureEdges.end()) {
			if ((iter->second).size() < circuitMinDegree) {

				//then this node will not contribute to the maximal mapping
				changeMade=true;
				
				//we will remove each edge one-by-one and as we remove, we will push our changes to the changes queue
				for (auto j=(iter->second).begin();j!=(iter->second).end();){
					architectureEdges[*j].erase(iter->first);
					changes.push(make_pair(*j,iter->first));
					j=(iter->second).erase(j);
				}
				iter = architectureEdges.erase(iter);
			} else{
				iter++;
			}
		}
	} while(changeMade);
	
	
	//returning all the edges we erased
	return make_pair(changes,architectureEdges);
}

 //rolling back changes if we went too far in pruning
 map<int,set<int>> rollBack(queue<pair<int,int>>changes, map<int,set<int>>architectureEdges) {
	 while (changes.size() > 0) {
		 pair<int, int> item = changes.front();
		 changes.pop();
		 int src = item.first;
		 int dst = item.second;
		 //re-adding edges
		 if (architectureEdges.find(src)==architectureEdges.end()){
			 set<int> edges;
			 edges.insert(dst);
			 architectureEdges[src]=edges;
		 } else{
			 architectureEdges[src].insert(dst);
		 }

		if (architectureEdges.find(dst)==architectureEdges.end()){
			 set<int> edges;
			 edges.insert(src);
			 architectureEdges[dst]=edges;
		 } else{
			 architectureEdges[dst].insert(src);
		 }
		 
		
	}
	return architectureEdges;

	//now all changes have been reverted
}

 //creating an initial array of distances of the shortest distance between any 2 physical qubits
 vector<vector<int>> getPhysicalDistancesFloydWarshall(map<int,set<int>> architectureEdges) {
		//we run the vector based floyd warshall algorithm here for shortest distances	
	 //initially setting max to be # of edges in total
	 int numEdges = 0;
	 for (auto i = architectureEdges.begin();i != architectureEdges.end();i++) {
		 numEdges += (i->second).size();
	 }
	 vector<vector<int>> distances(architectureEdges.size(), vector<int>(architectureEdges.size(),numEdges));

	 //setting up the algorithm
	 for (auto i = architectureEdges.begin();i != architectureEdges.end();i++) {
		 distances[i->first][i->first] = 0;
		 for (auto j = (i->second).begin();j != (i->second).end();j++) {
			 //distance 1 edges
			 distances[i->first][*j] = 1;
		 }
	 }

	 //running full algo here:
	 for (unsigned int i = 0;i < architectureEdges.size();i++) {
		 for (unsigned int j = 0; j < architectureEdges.size();j++) {
			 for (unsigned int z = 0; z < architectureEdges.size();z++) {
				 if (distances[j][z] > distances[j][i] + distances[i][z]) {
					 distances[j][z] = distances[j][i] + distances[i][z];
				 }
			}
		 }
	 }

	 //returning distance matrix
	 return distances;
}

 //naive distance implementation for now for distance between 2 physical qubits
 int getNaiveDistance(map<int, set<int>> architectureEdges, int Q1, int Q2) {
		//we basically just use bfs and get shortest distance
	 //queue<pair<int,int>> distanceQueue;
	 queue<int> distanceQueue;
	 queue<int> tempQueue;
	 set<int> closedSet;
	 int currDistance = 0;
	 //pushing source into queue
	 distanceQueue.push(Q1);
	 while (true) {
		 while (distanceQueue.size() > 0) {
			 int item = distanceQueue.front();
			 distanceQueue.pop();
			 if (closedSet.find(item) != closedSet.end()) {
				 //we already saw this
				 continue;
			 }

			 if (item == Q2) {
				 //we found the target
				 return currDistance;
			 }

			 //pushing all unseen neighbors to temp queue
			 for (auto i = architectureEdges[item].begin();i != architectureEdges[item].end();i++) {
				 if (closedSet.find(*i) == closedSet.end()) {
					 //then we can push this
					 tempQueue.push(*i);
				 }
			 }
		 }
		 //incrementing distance for next while loop
		 currDistance++;
		 //swapping queues now
		 queue<int> empty = distanceQueue;
		 distanceQueue = tempQueue;
		 tempQueue = empty;
	 }
	 
	 //for compilation
	 return 0;
 }

 //used alongside floyd warshall algorithm to generate the distance matrix
 int getSmartDistance(vector<vector<int>> distanceMatrix, int Q1, int Q2) {
	 return distanceMatrix[Q1][Q2];
 }

 int getDistance(vector<vector<int>> distanceMatrix, map<int,set<int>> architectureEdges, int Q1, int Q2) {
	 if (SMART) {
		 return getSmartDistance(distanceMatrix, Q1, Q2);
	 }
	 else {
		 return getNaiveDistance(architectureEdges, Q1, Q2);
	 }
 }



//helper to get my estimated partial mapping cost
 int partialMappingCost(map<int,int> mapping, set<GateNode*> remainingGates, vector<vector<int>> distanceMatrix, map<int,set<int>> architectureEdges) {
	//firstly pushing all gates to queue and addding them to closed set	
	 queue<GateNode*> remainingDependencies;
	 set<GateNode*> closedSet;
	 for (auto i = remainingGates.begin();i != remainingGates.end();i++) {
		 remainingDependencies.push(*i);
		 closedSet.insert(*i);
	 }

	 //now going gate by gate and inserting children as we go
	 int cost = 0;
	 while (remainingDependencies.size() > 0) {
		 //we will only consider cost for cnot gates
		 GateNode* toAppraise = remainingDependencies.front();
		 remainingDependencies.pop();
		 if (toAppraise->control == -1) {
			 continue;
		 }
		 if (mapping.find(toAppraise->control) != mapping.end() && mapping.find(toAppraise->target) != mapping.end()) {
			 //then both target and control are mapped to some physical qubit already
				//we increment cost by distance between their mapped physical qubits -1
				//notice that this means that if the qubits are adjacent, then cost is incremented by 0, which means the gate can get right
					//on its work with execution
			 cost += getDistance(distanceMatrix, architectureEdges ,mapping[toAppraise->control], mapping[toAppraise->target]) - 1;
		 }
		 else if (mapping.find(toAppraise->control) != mapping.end()) {
			 //then only control is mapped to some physical qubit
			 //we check if all the neighbors of control in the hardware mapping are mapped already
				//then we need to consider at least 1 swap in the cost (because there is no room to just assign the target mapping)
			 bool allMapped = true;
			 for (auto i = architectureEdges[toAppraise->control].begin();i != architectureEdges[toAppraise->control].end();i++) {
				 if (mapping.find(*i) == mapping.end()) {
					//then some neighbors is unmapped
					 allMapped = false;
					 break;
				}
			 }
			 if (allMapped) {
				 cost += 1;
			 }
		 }
		 else if (mapping.find(toAppraise->target) != mapping.end()) {
			 //then only target is mapped to some physical qubit
			 //we check if all neighbors of the target in the hardware mapping are mapped already
				//then we need to consider at least 1 swap in the cost (because there is no room to just assign the control mapping)
			 bool allMapped = true;
			 for (auto i = architectureEdges[toAppraise->target].begin();i != architectureEdges[toAppraise->target].end();i++) {
				 if (mapping.find(*i) == mapping.end()) {
					//then some neighbors is unmapped
					 allMapped = false;
					 break;
				}
			 }
			 if (allMapped) {
				 cost += 1;
			 }
		 }

		 if (toAppraise->controlChild != NULL ) {
				//then we need to check if target parent is in the closedSet or not for insertion (so we enforce level order)
				if (closedSet.find((toAppraise->controlChild)->targetParent) != closedSet.end() || remainingGates.find((toAppraise->controlChild)->targetParent)==remainingGates.end()) {
					//then we can add this child to the queue
					remainingDependencies.push(toAppraise->controlChild);
				}
			}

			if (toAppraise->targetChild != NULL && toAppraise->targetChild != toAppraise->controlChild) {
				//then we need to check if control parent is in the closedSet or not for insertion (so we enforce level order)
				if (closedSet.find((toAppraise->targetChild)->controlParent) != closedSet.end() || remainingGates.find((toAppraise->targetChild)->controlParent)!=remainingGates.end()) {
					//then we can add this child to the queue
					remainingDependencies.push(toAppraise->targetChild);
				}
			}

			closedSet.insert(toAppraise);
	 }


	 //returning my estimated cost
	 return cost;
 }

 //helper to get my estimated number of swaps between partial mappings
	//this will also help in my A* stitching algorithm
 int partialMappingStitchingCost(map<int,int> firstMap, map<int,int> secondMap,vector<vector<int>> distanceMatrix, map<int,set<int>> architectureEdges) {
	//all we do is get the maximum distance between any 2 mapped logical qubits in both mappings	
	 int maxDistance = 0;
	 auto iter = firstMap.begin();
	 while (iter != firstMap.end()) {
		 if (secondMap.find(iter->first) != secondMap.end()) {
			 //then this logical qubit (iter->first) is in both partial mappings
			 int calculatedDistance = getDistance(distanceMatrix, architectureEdges, firstMap[iter->first], secondMap[iter->first]);
			 if ( calculatedDistance> maxDistance) {
				 maxDistance = calculatedDistance;
			 }
		}
		 iter++;
	 }
	 //now we have the minimum # (lower bound) on swaps needed in order to transform a mapping into another mapping
	 return maxDistance;

 }

//returns only a single pathway for swaps to take instead of all of them
	//this is fine for the fast approach, but for the latency oriented approach we need to get all the swap pathways
 void singleEntryStitchMappings(pair<map<int, int>, set<int>> mapOne, pair<map<int, int>, set<int>> mapTwo, map<int, set<int>> architectureEdges) {
	//will implement later 
 }
/*
 pair<map<int, int>, vector<pair<int, int>>>improvedStitchMappings(pair<map<int, int>,set<int>> firstMap, pair<map<int, int>,set<int>> secondMap, map<int, set<int>> architectureEdges) {
	 map<int, int> mapOne = firstMap.first;
	 map<int, int> mapTwo = secondMap.first;

	 set<int> firstArch = firstMap.second;
	 set<int> secondArch = secondMap.second;
	//doing prep work first
		//we need to ensure that both sides of a swap pathway are mapped
		//otherwise we will not get an optimal number of swaps
		//consider LNN architecture and it may be seen that we might need to swap from either side to get a good solution
	 queue<map<int, int>> mappings;
	 bool iterationsNeeded = false;
	 for (auto z = mapTwo.begin();z != mapTwo.end();z++) {
		 if (mapOne.find(z->first) == mapOne.end()) {
			 iterationsNeeded = true;
			 break;
		 }
	 }

	 if (iterationsNeeded) {
		 mappings.push(mapOne);
		 //need to go through possible mappings for each unmapped qubit
		 for (auto i = 0;i != numLogical;i++) {
			 if (mapOne.find(i) == mapOne.end()) {
				 //need to map this

			 }
		 }
	 }

	 priority_queue<PriorityData> swapQueue;


	 //the big idea is that we consider 1 swap from every mapped qubit in the first mapping. we enqueue this and eventually we will hit mapping 2 with optimal swaps between these mapping due to cost function being an estimate
	 //initially just pushing the initial values with given cost 
	 int initialCost = partialMappingStitchingCost(mapOne, mapTwo, distanceMatrix, architectureEdges);
	 swapQueue.push(PriorityData(initialCost, vector<pair<int, int>>(), mapOne, mapOne));

	 while (swapQueue.size() > 0) {
		 PriorityData current = swapQueue.top();
		 swapQueue.pop();
		 map<int, int> transformedMapping = current.transformedMapping;
		 map<int, int> initialMapping = current.initialMapping;
		 vector<pair<int, int>> swaps = current.swapsTaken;
		 //checking if stuff is satisfied for mapping
		 for (auto z = mapTwo.begin();z != mapTwo.end();z++) {
			 if (transformedMapping.find(z->first) != transformedMapping.end() && transformedMapping[z->first] == mapTwo[z->first]) {
				 //then we found our swap pathway!
				 return make_pair(initialMapping, swaps);
			 }
		 }

		 //going through every mapped qubit and considering exactly 1 swap
		 for (auto )
	 }
	//if a qubit is mapped in first mapping, but not in second, no worries
	//if a qubit is mapped in second, but not in first, need to go through possibilities of what qubit could be beforehand and enqueue
	 //set of unmapped physical qubits that are neighbors to currently mapped physical qubits
	 set<int> unmappedNeighbors;
	 //logical qubits that need some initial mapping right now
	 vector<int> needMapping;
	 for (auto i = secondMap.begin();i != secondMap.end();i++) {
		 if (firstMap.find(i->first) == firstMap.end()) {

		 }
	}

	 while (swapQueue.size() > 0) {
		 //pop lowest estimated cost swap pathway
		 //consider a single swap and then 
	}

	 return pair<map<int, int>, vector<pair<int, int>>>();
 }
 */
 

 //function to get a list of swaps needed to transform one mapping to another
	//returns as the first element of the pair the modified initial mapping (additional logical qubits mapped to mapOne, if necessary)
	//second element of the returned pair is the list of swaps to take
 pair<map<int,int>,vector<pair<int,int>>> stitchMappings(pair<map<int,int>,set<int>> mapOne, pair<map<int,int>,set<int>> mapTwo, map<int,set<int>> architectureEdges) {
	//initializing a priority queue
	
	 //priority queue of the next swap to take
	 //the queue is the queue (order) of swaps to take
		//the map is the modified mapping with the swap
	 //queue structure of tuple	
		//1) this is the cost function for pursuing this mapping, (cost of swaps taken + heuristic cost on transformed mapping)
		//2) this is the initial mapping, this is provided in case changes need to be made 
			//(think about a case where the initial mapping has an unmapped logical qubit where the target mapping has mapped that logical qubit)
		//3) this is the queue of swaps taken (logicalQ_1, logicalQ_2)
		//4) this is the transformed mapping 
	 priority_queue<PriorityData> swapQueue;
	 /*
	 priority_queue<tuple<int,
		 vector<pair<int,int>>,
		 map<int,int>,
		 map<int,int>>,vector<tuple<int,
		 vector<pair<int,int>>,
		 map<int,int>,
		 map<int, int>>>, > swapQueue;
		 */
	 map<int, int> firstMap = mapOne.first;
	 map<int, int> transformedMap = firstMap;
	 set<int> firstArchMap = mapOne.second;
	 map<int, int> secondMap = mapTwo.first;
	 set<int> secondArchMap = mapTwo.second;

	 for (auto j = transformedMap.begin();j != transformedMap.end();j++) {
		 //printf("INITIAL MAPPING %d mapped to %d\n", j->first, j->second);
	 }
	

	 //enqueueing the base state
	 //printf("ENQUEING BASE STATE!\n");
	 swapQueue.push(PriorityData(partialMappingStitchingCost(firstMap, secondMap,distanceMatrix,architectureEdges), vector<pair<int, int>>(), firstMap, transformedMap));
	 int iterCount = 0;
	 while (swapQueue.size() > 0) {
		 //printf("SIZE OF SWAP QUEUE: %d\n", swapQueue.size());
		 //printf("GOT INSIDE!\n");
		//getting the lowest cost path from the queue
		 PriorityData result = swapQueue.top();
		 //tuple<int,vector<pair<int, int>>,map<int, int>,map<int, int>>  result = swapQueue.top();
		 swapQueue.pop();
		 //vector<pair<int, int>> swapsTaken = get<1>(result);
		 //map<int, int> initialMapping = get<2>(result);
		 //map<int, int> transformedMapping = get<3>(result);
		 //int currCost = get<0>(result);
		 vector<pair<int, int>> swapsTaken = result.swapsTaken;
		 map<int, int> initialMapping = result.initialMapping;
		 map<int, int> transformedMapping = result.transformedMapping;
	     int currCost = result.cost;		 
		 //printf("CURRENT TRANSFORMED MAPPING:\n");
		 for (auto o = transformedMapping.begin();o != transformedMapping.end();o++) {
			 //printf("%d mapped to %d, ", o->first, o->second);
		}
		 if (swapsTaken.size() != 0) {
			 //printf("LAST CHOSEN SWAP WAS BETWEEN %d and %d with Cost %d\n", swapsTaken[swapsTaken.size() - 1].first, swapsTaken[swapsTaken.size() - 1].second, currCost);
		 }
		 //printf("WITH COST: %d", currCost);
		 //printf("\n");
		 //debug
		 
		 //checking equality 
		 bool equal = true;
		 for (auto h = secondMap.begin();h != secondMap.end();h++) {
			 if (transformedMapping.find(h->first) == transformedMapping.end()) {
				 //printf("HIT THIS ERROR!\n");
				 //need to check if we can just assign this to all such qubits. if so, then we are done
					//i.e, if everything is mapped to everything it is supposed to, except 2 qubits here, then we can just map them and we are done
				 bool allElseMapped = true;
				 for (auto z = secondMap.begin();z != secondMap.end();z++) {
					 if (transformedMapping.find(z->first) != transformedMapping.end()) {
						 if (transformedMapping[z->first] != z->second) {
							 //printf("LOGICAL %d mapped to %d in trasnformed and %d in second\n", z->first, transformedMapping[z->first], z->second);
							 allElseMapped = false;
							 break;
						 }
					 }
				 }

				 if (allElseMapped) {
					 //printf("ALL ELSE MAPPED!\n");
					 //we will try to artificially create the mapping here because it may be possible
					 map<int, int> newInitialMapping = initialMapping;
					 map<int, int> newTransformed = transformedMapping;
					 bool allMapped = true;
					 for (auto z = secondMap.begin();z != secondMap.end();z++) {
						 if (transformedMapping.find(z->first) == transformedMapping.end()) {
							 //lets see if we can map this
							 for (auto b = transformedMapping.begin();b != transformedMapping.end();b++) {
								 if (b->second == z->second) {
									 //not possible
									 allMapped = false;
									 break;
								 }
								 else {
									 //lets add this to our initial mapping and update transformed
									 newInitialMapping[z->first] = z->second;
									 newTransformed[z->first] = z->second;
								 }
							 }
						 }
					 }
					 if (allMapped) {
						 //mapping was possible!
						 return make_pair(newInitialMapping, swapsTaken);
					 }
				 }

				 
				 //then automatically not equal
				 equal = false;
				 break;
			 }
			 if (transformedMapping[h->first] != h->second) {
				 //mapping not satisfied
				 equal = false;
				 break;
			 }
		 }
		 if (equal) {
			 //then we hit the end
			 //returning the initial mapping, and list of swaps to take
			 //printf("FOUND MAPPING PATHWAY!\n");
			 return make_pair(initialMapping, swapsTaken);
		 }

		 //else we need to think in terms of 1 swap granularities
		 //we have some cases to think about:
			//1) suppose a qubit exists in mapping 1 but not in mapping 2
				//this is no issue because then this qubit does not contribute to the next section of the dependence graph
			//2) suppose a logical qubit is mapped in mapping 2 but not in mapping 1
				//then that logical qubit must be involved in an adjacent cnot gate
				//therefore, we try mapping common qubits first, and as we consider swaps with unmapped neighbors, we assign qubits
			//3) if a logical qubit is mapped in 1 and 2
				//then we just consider possible swaps of distance 1 away
		 for (auto j = transformedMapping.begin();j != transformedMapping.end();j++) {
			 if (secondMap.find(j->first) == secondMap.end()) {
				 //then this qubit does not contribute to the next maximal mapping 
				 //printf("LOGICAL QUBIT %d DOES NOT CONTRIBUTE TO NEXT MAPPING!\n", j->first);
				 continue;
			 }
			 if (transformedMapping[j->first] != secondMap[j->first]) {
				 //printf("CASE 1 NEED TO SWAP!\n");
				 //will need to consider 1 swap applied to j->first or 1 swap applied to the logical qubit mapped to secondMap[j->first]
				 // for loop below considering a single swap applied to logical qubit mapping j->first
				 for (auto i = architectureEdges[transformedMapping[j->first]].begin();i != architectureEdges[transformedMapping[j->first]].end();i++) {
					 //need to check if this qubit is mapped to in the current mapping or not
						//if so, we consider a swap with the associated logical qubit
						//if not, then we check if there are qubits that are mapped in the result that are not mapped in current
							//this is a possible assignment for them
					 bool nodePresent = false;
					 int associatedLogical = -1;
					 for (auto z = transformedMapping.begin();z != transformedMapping.end();z++) {
						 //checking if *i is in here
						 if (z->second == *i) {
							 nodePresent = true;
							 associatedLogical = z->first;
							 break;
						 }
					 }

					 if (nodePresent) {
						 //printf("REACHED EEEE!\n");
						//then the edge is mapped in the current mapping to logical qubit z
						//we modify our swaps, get the transformed mapping, update the cost and push back to our priority queue
						 vector<pair<int, int>> newSwaps = swapsTaken;
						 if (newSwaps.size() > 0) {
							
							pair<int, int> lastSwapTaken = newSwaps[newSwaps.size() - 1];
							 //checking  if the swap is redundant (reverses the previous swap)
							// then this swap is useless and we can continue
							if ((lastSwapTaken.first == associatedLogical && lastSwapTaken.second == j->first) || (lastSwapTaken.first == j->first && lastSwapTaken.second == associatedLogical)) {
								//then the swap is redundant
								continue;
							}
						 }
						 //adding proposed swap to vector of swaps
						 newSwaps.push_back(make_pair(associatedLogical, j->first));
						 //printf("CONSIDERING SWAP BETWEEN %d and %d\n", associatedLogical, j->first);
						
						 //getting the transformed mapping
						 map<int, int> newTransformedMapping = transformedMapping;

						 //swapping the logical qubit assignments
						 newTransformedMapping[associatedLogical] = transformedMapping[j->first];
						 newTransformedMapping[j->first] = *i;

						 //debug
						 //printf("CHECKING TRANSFORMED MAPPING!\n");
						 for (auto p = newTransformedMapping.begin();p != newTransformedMapping.end();p++) {
							 //printf("%d mapped to %d, ", p->first, p->second);
						 }
						 //printf("\n");

						 //getting new cost of the swap pathway
						 int newCost = newSwaps.size() + partialMappingStitchingCost(newTransformedMapping, secondMap, distanceMatrix, architectureEdges);
						 //printf("COST ASSIGNED OF %d\n", newCost);
						 //pushing to queue
						 swapQueue.push(PriorityData(newCost, newSwaps, initialMapping, newTransformedMapping));
					 }
					 else {
						 //printf("CONSIDER SWAP FOR AN UNMAPPED NODE CASE 1\n");
						//then the considered swap has a node that is not considered in the current mapping
							//if there is a logical in the target mapping that is not mapped in the current mapping, we can modify our initial mapping here
						      //otherwise, we continue
						 for (auto b = secondMap.begin();b != secondMap.end();b++) {
							 if (transformedMapping.find(b->first) == transformedMapping.end()) {
								 //then this is a logical qubit involved in the target mapping that is not in the current mapping
								 //we have to modify our initial mapping to account for this
								 map<int, int> newInitialMapping = initialMapping;
								 newInitialMapping[b->first] = *i;
								 //printf("CHANGED INITIAL MAPPING -1\n");
								 //updating swaps
								 vector<pair<int, int>>  newSwaps = swapsTaken;
								 newSwaps.push_back(make_pair(b->first, j->first));
								 //printf("CONSIDERING SWAP BETWEEN %d and %d\n", b->first, j->first);

								 //transforming mapping
								 map<int, int> newTransformedMapping = transformedMapping;
								 newTransformedMapping[b->first] = transformedMapping[j->first];
								 newTransformedMapping[j->first] = *i;

								 //updating cost
								 int newCost = newSwaps.size() + partialMappingStitchingCost(newTransformedMapping, secondMap, distanceMatrix, architectureEdges);
								 //printf("WITH COST %d\n", newCost);
								 //pushing to queue
								 swapQueue.push(PriorityData(newCost, newSwaps, newInitialMapping, newTransformedMapping));
							}
						}

					 }
				 }
				 
				 //checking secondMap[j->first]
					//specifically, is there a logical qubit in transformedMapping that is mapped to the same physical qubit as j->first in the target mapping?
					//if so, then the swap to consider is easy
					//if not, then we might need to change our initial mapping to account for this

				 //for loop to check if there is a logical qubit in the first mapping that is mapped to the same physical qubit as j->first in the second mapping
				 int associatedLogicalQubit = -1;
				 for (auto i = transformedMapping.begin(); i != transformedMapping.end();i++) {
					 if (i->second == secondMap[j->first]) {
						 //then we found it
						 associatedLogicalQubit = i->first;
						 break;
					 }
				 }
				
				 
				 if (associatedLogicalQubit != -1) {
					 //printf("GOT HERE! with associatedLogicalQubit %d\n",associatedLogicalQubit);
					 //then the physical qubit mapped to j->first in the target mapping has an associated logical qubit in the current mapping
						//we can just consider all swaps on associatedLogicalQubit of distance 1 now
					 for (auto z = architectureEdges[transformedMapping[associatedLogicalQubit]].begin();z != architectureEdges[transformedMapping[associatedLogicalQubit]].end();z++) {	
						 int receivingSwapQubit = -1;
						 for (auto c = transformedMapping.begin();c != transformedMapping.end();c++) {
							if (c->second == *z) {
								receivingSwapQubit = c->first;
								break;
							}
						 }
						 if (receivingSwapQubit==-1) {
							 //printf("OOOOH HERE!\n");
							 //then the receiving physical qubit is not mapped by some logical qubit in the transformed Mapping 
							 //we need to check if it needs to be assigned a mapping
							 for (auto b = secondMap.begin();b != secondMap.end();b++) {
								 if (transformedMapping.find(b->first) == transformedMapping.end()) {
									 //then this is a logical qubit involved in the target mapping that is not in the current mapping
									 //we have to modify our initial mapping to account for this
									 //printf("INITIAL MAPPING CHANGED 1!\n");
									 map<int, int> newInitialMapping = initialMapping;
									 newInitialMapping[b->first] = *z;
									 //updating swaps
									 vector<pair<int, int>>  newSwaps = swapsTaken;
									 newSwaps.push_back(make_pair(associatedLogicalQubit, b->first));
									 //printf("CONSIDERING SWAP BETWEEN %d and %d\n", associatedLogicalQubit, b->first);

									 //transforming mapping
									 map<int, int> newTransformedMapping = transformedMapping;
									 newTransformedMapping[b->first] = transformedMapping[associatedLogicalQubit];
									 newTransformedMapping[associatedLogicalQubit] = b->second;

									 //updating cost
									 int newCost = newSwaps.size() + partialMappingStitchingCost(newTransformedMapping, secondMap, distanceMatrix, architectureEdges);
									 //printf("WITH COST %d\n", newCost);
									 //pushing to queue
									 swapQueue.push(PriorityData(newCost, newSwaps, newInitialMapping, newTransformedMapping));
								}
							}	 
						 } else {
							 //printf("OVER HEREEE!\n");
							 //then we can just consider a swap between associatedLogicalQubit and receivingSwapQubit
								//we should check first if the chosen swap was already done previously (dont want to be redundant)
							 vector<pair<int, int>> newSwaps = swapsTaken;
							 if (newSwaps.size() > 0 && (
								 (newSwaps[newSwaps.size() - 1].first == associatedLogicalQubit && newSwaps[newSwaps.size() - 1].second == receivingSwapQubit) ||
								 (newSwaps[newSwaps.size() - 1].first == receivingSwapQubit && newSwaps[newSwaps.size() - 1].second == associatedLogicalQubit))) {
							 
								 //then this is a redundant swap
								 continue;
							 }
							 newSwaps.push_back(make_pair(associatedLogicalQubit, receivingSwapQubit));
							 //printf("CONSIDERING SWAP 1 BETWEEN %d and %d\n", associatedLogicalQubit, receivingSwapQubit);

							 //transforming the mapping
							 map<int, int> newTransformedMapping = transformedMapping;
							 newTransformedMapping[associatedLogicalQubit] = transformedMapping[receivingSwapQubit];
							 newTransformedMapping[receivingSwapQubit] = transformedMapping[associatedLogicalQubit];

							 //getting new cost
							 int newCost = newSwaps.size() + partialMappingStitchingCost(newTransformedMapping, secondMap, distanceMatrix, architectureEdges);
							 //printf("WITH COST %d\n", newCost);

							 //pushing to priority queue
							 swapQueue.push(PriorityData(newCost, newSwaps, initialMapping, newTransformedMapping));
						 }
					 }
				 }
				 else {
					 //printf("OH NO!\n");
					 //then we are in the case where target mapping has mapped a physical qubit that is not mapped to in the current mapping
						//we might need to change the initial mapping here, if possible
						//we might need to change the initial mapping here twice (for considering swaps of distance 1 where the target is not mapped!
					 for (auto z = secondMap.begin();z != secondMap.end();z++) {
						 if (transformedMapping.find(z->first) == transformedMapping.end()) {
							 //we found a logical qubit in the second mapping that is not mapped in the current mapping
								//we can assign it the physical qubit  secondMap[j->first] in the initial mapping
							 //printf("INITIAL MAPPING CHANGED 2\n");
							 map<int, int> newInitialMapping = initialMapping;
							 newInitialMapping[z->first] = secondMap[j->first];
							 map<int, int> newTransformedMapping = transformedMapping;
							 newTransformedMapping[z->first] = secondMap[j->first];
							 //now considering all swaps of 1 distance from this logical qubit
							 for (auto c = architectureEdges[secondMap[j->first]].begin();c != architectureEdges[secondMap[j->first]].end();c++) {
								// seeing if the target edge is mapped by any logical qubit in transformed mapping
									//if so, we consider a swap between logical qubits
									//if not, then we need to modify our initial mapping
								 int associatedQubit = -1;
								 for (auto d = transformedMapping.begin();d != transformedMapping.end();d++) {
									 if (d->second == *c) {
										 //then we found the logical qubit associated to the above physical qubit
										 associatedQubit = d->first;
										 break;
									 }
								}
								
								 if (associatedQubit == -1) {
									//then we modify our initial mapping further if possible
									 for (auto d = secondMap.begin();d != secondMap.end();d++) {
										 if (newTransformedMapping.find(d->first) == newTransformedMapping.end()) {
											 //then we found a possible initial mapping modifier
											 newInitialMapping[d->first] = *c;
											
											 //considering the swap
											 vector<pair<int, int>> newSwaps = swapsTaken;
											 newSwaps.push_back(make_pair(d->first, z->first));
											 //printf("CONSIDERING SWAP 2 BETWEEN %d and %d\n", d->first, z->first);

											 newTransformedMapping[z->first] = *c;
											 newTransformedMapping[d->first] = secondMap[j->first];

											 //getting new cost
											 int newCost = newSwaps.size() + partialMappingStitchingCost(newTransformedMapping, secondMap, distanceMatrix, architectureEdges);
											 //printf("WITH COST %d\n", newCost);

											 swapQueue.push(PriorityData(newCost, newSwaps, newInitialMapping, newTransformedMapping));
										 }
									}
								 }
								 else {
									 //now we just consider a swap between these qubits
									 vector<pair<int, int>> newSwaps = swapsTaken;
									 newSwaps.push_back(make_pair(associatedQubit, z->first));
									 //printf("CONSIDERING SWAP 3 BETWEEN %d and %d\n", associatedQubit, z->first);

									 //updating our transformed mapping
									 newTransformedMapping[associatedQubit] = secondMap[j->first];
									 newTransformedMapping[z->first] = transformedMapping[associatedQubit];

									 //getting new cost
									 int newCost = newSwaps.size() + partialMappingStitchingCost(newTransformedMapping, secondMap, distanceMatrix, architectureEdges);
									 //printf("WITH COST %d\n", newCost);
									 //pushing this search node into the queue
									 swapQueue.push(PriorityData(newCost, newSwaps, newInitialMapping, newTransformedMapping));
								 }
							 }

						 }
					 }
				 }
				  
			 }
		 }
		 iterCount++;
	 }
	 return make_pair(map<int,int>(), vector<pair<int,int>>());
 }

 //goes through a pathway of swaps and returns a list of sets, where each set contains all the swaps that can occur at once
	//this will help us generate circuits with lower hits to latency due to swaps
 vector<set<pair<int,int>>> parallelizeSwapPathway(vector<pair<int,int>> swapPathway) {
	//idea: for each swap, go through every other swa	
	 //idea: use a set to determine conflicts (since these swaps are in order)
			//basic loop idea: go through the entire pathway, insert any swaps that do not conflict into a set
	 vector<set<pair<int, int>>> parallelizedSwaps;
	 set<pair<int, int>> addingSet;
	 set<int> currentSet;
	 while (true) {
		 if (swapPathway.size() == 0) {
			 //we are done
			 break;
		 }
		 for (auto i = swapPathway.begin();i!=swapPathway.end();) {
			 pair<int, int> swap = *i;
			 if (currentSet.find(swap.first) == currentSet.end() && currentSet.find(swap.second) == currentSet.end()) {
				 //then this swap does not conflict with the previous ones
				 addingSet.insert(swap);
				 currentSet.insert(swap.first);
				 currentSet.insert(swap.second);
				 i = swapPathway.erase(i);
			 }
			 else {
				 i++;
			 }
		 }
		 parallelizedSwaps.push_back(addingSet);
		 //resetting the set
		 addingSet = set<pair<int, int>>();
		 //resetting the set of qubits used
		 currentSet = set<int>();

		 
	 }

	 //now parallelizedSwaps has a vector where each set is a collection of swaps that can be executed in parallel
	 return parallelizedSwaps;
 }

 //helper function that , when given a list of starting gates, and a list of ending gates, returns a queue of gates up to the ending gates point
 queue<GateNode*> getSatisfiedGates(set<GateNode*> beginGates, set<GateNode*> endGates) {
	 set<GateNode*> start = beginGates;
	 set<GateNode*> closedSet;
	 queue<GateNode*> orderQueue;
	 queue<GateNode*> finalQueue;

	 while (start.size() > 0) {
			////printf("current size: %d\n",start.size());
			for (auto i = start.begin();i!=start.end();) {
				
				if (start.find((*i)->controlParent) == start.end() && start.find((*i)->targetParent) == start.end()) {
					//then we add it into queue to get a level order
					orderQueue.push((*i));
					i=start.erase(i);
					////printf("gate pushed!\n");	
				} else{
					i++;
				}
			}
	 }
	 
	 while (orderQueue.size() > 0) {
		 GateNode* toConsider = orderQueue.front();
		 orderQueue.pop();
		 if (endGates.find(toConsider) != endGates.end()) {
			 //this gate was not satisfied
			 continue;
		 }

		 if (closedSet.find(toConsider) != closedSet.end()) {
			 //we already saw this gate
			 continue;
		 }

		 //otherwise we add children to order queue and we can add this gate to the print queue
		 finalQueue.push(toConsider);
		 closedSet.insert(toConsider);
		  
		 if (toConsider->controlChild!=NULL){
					GateNode* child = toConsider->controlChild;
					if (child->control == toConsider->control){
						//then we need to check target parent
						if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end() || start.find(child->targetParent)==start.end()){
							orderQueue.push(child);
						}
					} else {
						//then we need to check control parent
						if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end() || start.find(child->controlParent)==start.end()){
							orderQueue.push(child);
						}
					}
				}

				if (toConsider->targetChild!=NULL && toConsider->targetChild != toConsider->controlChild){
					GateNode* child = toConsider->targetChild;
					
					if (child->control == toConsider->target){
						//then we need to check target parent
						if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end() || start.find(child->targetParent)==start.end()){
							orderQueue.push(child);
						}
					} else {
						
						//then we need to check control parent
						if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end() || start.find(child->controlParent)==start.end()){
							orderQueue.push(child);
						}
					}
				}
	 }

	 //returning the final print queue
	 return finalQueue;
 }

 //iterates through every possibility of unmapped   
 /*
 tuple<map<int, int>, set<int>> getBestResult(tuple<map<int, int>, set<int>, set<GateNode*>>) {

 }
 */
//lazy function to insert swaps (just tacks on parallelized swaps at the end of every maximal partition)
 void betterLazySwapCircuitBuilder(map<int,set<int>> architectureEdges, set<GateNode*> startGates, vector<vector<int>> distanceMatrix) {
 //since our initial mapping varies until everything is mapped, we have to print out gates at the end
 //but the process is this: find the best initial maximal mapping (maximal Mapping cost)
	 //gates satisfied for every step of the maximal mapping
	 vector<queue<GateNode*>> maximalMappingSatisfiedGates;
	 //swaps that need to be taken after every maximal partition, if any
	 vector<vector<pair<int, int>>> levelSwaps;
	 //maintaining our initial mapping (every qubit may not be mapped at the start)
	 map<int, int> initialMapping;
	 //lastMapping is tracked in case we need to perform swaps
	 pair<map<int, int>,set<int>> lastMapping;
	 set<GateNode*> currentGates = startGates;
	 while (currentGates.size() > 0) {
		 //doing maximal mapping first
		 vector<tuple<map<int,int>,set<int>,set<GateNode*>>> partialMapping=maximalMapper(architectureEdges, currentGates);

		

		 
		 //if this is first maximal mapping we will need to find the best one for remaining cost and set that as our initial mapping
		 if (lastMapping.first.size() == 0) {
			 //finding the best mapping out of these
			 set<GateNode*> endGates;
			 int lowestCost = -1;
			 for (auto z = partialMapping.begin();z != partialMapping.end();z++) {
				 int cost = partialMappingCost(get<0>(*z), get<2>(*z), distanceMatrix, architectureEdges);
				if (lowestCost == -1 || cost < lowestCost) {
					lastMapping = make_pair(get<0>(*z), get<1>(*z));
					lowestCost = cost;
					endGates = get<2>(*z);
				} 
			 }

			 //setting our initial mapping to last mapping
			 initialMapping = lastMapping.first;

			 //adding gates to the print queue
			 queue<GateNode*> toPrint = getSatisfiedGates(currentGates, endGates);

			 //updating currentGates
			 currentGates = endGates;

			//pushing print queue to vector
			maximalMappingSatisfiedGates.push_back(toPrint);
			printf("FINISHED FIRST ITERATION!\n");
			continue;
		 }
		 else {
			 printf("STARTED OTHER ITERATION!\n");
			//need to find the closest mapping to the previous mapping, and then consider swaps
			 pair<map<int, int>, set<int>> bestMapping;
			 int bestMappingCost = -1;
			 set<GateNode*> endGates;
			 for (auto z = partialMapping.begin();z != partialMapping.end();z++) {
				 int cost = partialMappingStitchingCost(lastMapping.first, get<0>(*z), distanceMatrix, architectureEdges);
				 if (bestMappingCost==-1 || cost<bestMappingCost) {
					 //then we can just set this as the best
					 bestMapping = make_pair(get<0>(*z), get<1>(*z));
					 bestMappingCost = cost;
					 endGates = get<2>(*z);
				}
			 }

			 //now we find the list of swaps to take
			 pair<map<int, int>, vector<pair<int, int>>> swaps = stitchMappings(lastMapping, bestMapping, architectureEdges);

			 //push back list of swaps to vector
			 levelSwaps.push_back(swaps.second);
			 
			 map<int, int> toModify = lastMapping.first;
			 //we update our initial mapping here, if necessary
			 for (auto z = swaps.first.begin(); z != swaps.first.end();z++) {
				 if (initialMapping.find(z->first) == initialMapping.end()) {
					 initialMapping[z->first] = z->second;
					 toModify[z->first] = z->second;
				 }
			 }

			 //updating our mapping with the swaps
			 for (auto z = swaps.second.begin();z != swaps.second.end();z++) {
				//swaps are between logical qubits
				 int tmp = toModify[z->first];
				 toModify[z->first] = toModify[z->second];
				 toModify[z->second] = tmp;
			 }

			 lastMapping = make_pair(toModify, set<int>());

			 //adding gates to the print queue
			 queue<GateNode*> toPrint = getSatisfiedGates(currentGates, endGates);

			 //updating currentGates
			 currentGates = endGates;

			//pushing print queue to vector
			maximalMappingSatisfiedGates.push_back(toPrint);
			printf("FINISHED NEXT ITERATION\n");
		 }



	 }

	 printf("REACHED HERE!\n");

	 //if there are unmapped qubits in the initial mapping, then they do not matter and we can assign them whatever we want here
		//specifically, we can assign them any unmapped physical qubit
	 while (initialMapping.size() != numLogical) {
		 for (int i = 0;i < numLogical;i++) {
			 if (initialMapping.find(i) == initialMapping.end()) {
				 //then we can assign them any unmapped physical qubit
				 for (auto j = architectureEdges.begin();j != architectureEdges.end();j++) {
					 bool mapped = false;
					 for (auto z = initialMapping.begin();z != initialMapping.end();z++) {
						 if (z->second == j->first) {
							 mapped = true;
							 break;
						}
					 }
					 if (!mapped) {
						 //then we can map it
						 initialMapping[i] = j->first;
					 }
				 }
				 break;
			 }
		 }
	 }
	 
	 //keeping a copy to print at the very end
	 map<int, int> copyInitialMapping = initialMapping;
	 //proceeding with printing/writing the circuit
		//we just start with the initial mapping --> print the first set of gates ---> apply swaps ---> print next set of gates, and so on
	 for (int i = 0;i != maximalMappingSatisfiedGates.size();i++) {
		 if (i != 0) {
			 //need to apply swaps first before printing
			 //we parallelize swaps and then print
			 //we have exactly one less swap entry than levels (start level does not need swaps)
			 vector<pair<int, int>> swaps = levelSwaps[i - 1];
			 vector<set<pair<int, int>>> parallel = parallelizeSwapPathway(swaps);
			 for (auto p = parallel.begin();p != parallel.end();p++) {
				 for (auto z = (*p).begin();z != (*p).end();z++) {
					 //adding a swap to initial mapping and also printing the swap
					 cout << "swp q[" << initialMapping[z->first] << "], q[" << initialMapping[z->second] << "]; // swap logical qubits q" << z->first << ", and q" << z->second << "\n";
					 //applying the swap
					 int tmp = initialMapping[z->first];
					 initialMapping[z->first] = initialMapping[z->second];
					 initialMapping[z->second] = tmp;
				 }
			 }
		 }
		 queue<GateNode*> satisfied = maximalMappingSatisfiedGates[i];
		 while (satisfied.size() > 0) {
			 GateNode* toPrint = satisfied.front();
			 satisfied.pop();
			 
			 //printing name and qubits involved
			 if (toPrint->control == -1) {
				 //printing single qubit gate
				 cout << toPrint->name << " q[" << initialMapping[toPrint->target] << "]; // formerly on qubit q[" << toPrint->target << "]\n";
			 }
			 else {
				 //printing 2 qubit gate
				 cout << toPrint->name << " q[" << initialMapping[toPrint->control] << "], q[" << initialMapping[toPrint->target] << "]; // formerly on qubits q[" << toPrint->control << "], q[" << toPrint->target << "]\n";
			 }
		 }
	 }

	 //printing the initial mapping
	 cout << "INITIAL MAPPING!: \n";
	 for (auto z = copyInitialMapping.begin();z != copyInitialMapping.end();z++) {
		 cout << "Logical qubit q" << z->first << " mapped to physical qubit Q" << z->second <<"\n";
	 }
 }

 //lazy function to insert swaps (just tacks them on at the end of every maximal partition)
 void lazySwapCircuitBuilder(map<int,set<int>> architectureEdges, set<GateNode*> startGates,vector<vector<int>> distanceMatrix) {
	 set<GateNode*> nextGates = startGates;
	 //set of gates to build on for building the actual circuit
	 //set<GateNode*> resultGates;
	 //managing our initial mapping
	 map<int, int> initialMapping;
	 pair<map<int, int>,set<int>> lastMapping;
	 int count = 0;
	 while (nextGates.size() > 0) {
		 //we start by getting a maximal mapping
		 vector<tuple<map<int, int>, set<int>, set<GateNode*>>> result = maximalMapper(architectureEdges, nextGates);
		 
		 //we go through all the mappings and pick the best partial mapping according to my partial mapping cost function and difference function
			//we first check the difference between the previous mapping, if not null
			//then the tie-breaker is partial mapping cost
			//if still ties, we just pick the first one
		 tuple<map<int, int>, set<int>, set<GateNode*>>  bestTuple;
		 if (get<0>(lastMapping).size() == 0) {
			 int currTupleCost = -1;
			 for (auto i = result.begin();i != result.end();i++) {
				 //going through entire vector and seeing which swap is the best
				 int cost = partialMappingCost(get<0>(*i), get<2>(*i), distanceMatrix, architectureEdges);
				 if (get<0>(bestTuple).size() == 0) {
					 bestTuple = *i;
					 currTupleCost = cost;
					 continue;
				 }
				 else if (cost < currTupleCost) {
					 //then we proceed with this mapping instead
					 bestTuple = *i;
					 currTupleCost = cost;
				 }
			 }
		 }
		 else {
			 //we will have to consider "closeness" to the last mapping along with predicted cost for the remaining circuit
				//so we do 2 passes in result: 1 for closeness to last mapping first
											  //2nd pass for predicted remaining cost if size is large

			 //doing first pass for closeness
			 vector<int> closenessValues;
			 vector<tuple<map<int, int>, set<int>, set<GateNode*>>> closestResults;
			 int lowestCloseness = -1;
			 for (auto z = result.begin();z != result.end();z++) {
				 int closenessResult = partialMappingStitchingCost(lastMapping.first, get<0>(*z), distanceMatrix, architectureEdges);
				 if (lowestCloseness == -1) {
					 lowestCloseness = closenessResult;
				 }
				 else if (closenessResult < lowestCloseness) {
					 lowestCloseness = closenessResult;
				 }
				 closenessValues.push_back(closenessResult);
			 }

			 //putting closest values into another vector
			 for (unsigned int z = 0;z != closenessValues.size();z++) {
				 if (closenessValues[z] == lowestCloseness) {
					 //then we know the associated positional mapping is at least "close enough" to our last mapping
					 closestResults.push_back(result[z]);
				 }
			 }

			 //now we have the closest values in another vector 
			 if (closestResults.size() == 1) {
				 //we found our match, we can set and break
				 bestTuple = closestResults[0];

			 }
			 else {
				 //need to do second pass but with partial mapping costs now, and we only pick one
				 int lowestMappingCost = -1;
				 for (auto k = closestResults.begin();k != closestResults.end();k++) {
					 int estimatedCostRemaining = partialMappingCost(get<0>(*k), get<2>(*k), distanceMatrix, architectureEdges);
					 if (lowestMappingCost == -1) {
						 bestTuple = *k;
						 lowestMappingCost = estimatedCostRemaining;
					 }
					 else if (estimatedCostRemaining < lowestMappingCost) {
						 bestTuple = *k;
						 lowestMappingCost = estimatedCostRemaining;
					 }
				 }
			 }



		 }
		 /*
		 if (count == 0) {
			 pair<map<int,int>, set<int>> result = getBestResult(bestTuple);
		 }
		 */
		 //getting swaps, if necessary
		 if (get<0>(lastMapping).size() == 0) {
			 //then there is no mapping to stitch to previously, we can replace
			 lastMapping = make_pair(get<0>(bestTuple), get<1>(bestTuple));
			 initialMapping = get<0>(bestTuple);
		 }
		 else {

		 
			 pair<map<int, int>, set<int>> currentMaximalMapping = make_pair(get<0>(bestTuple), get<1>(bestTuple));
			 //then we need to stich swaps to the previous mapping and adjust our last mapping value and initial mapping if necessary
			 pair < map<int, int>, vector<pair<int, int>>> stitchingResult = stitchMappings(lastMapping, currentMaximalMapping, architectureEdges);
			 //adjusting the initial mapping
			 map<int, int> newInitialMapping = stitchingResult.first;
			 map<int, int> newLastMapping = lastMapping.first;
			 for (auto j = newInitialMapping.begin();j != newInitialMapping.end();j++) {
				 //adjusting the initial mapping as we map more qubits
				 if (initialMapping.find(j->first) == initialMapping.end()) {
					 //we should add this to our initial mapping and newLastMapping
					 initialMapping[j->first] = j->second;
					 newLastMapping[j->first] = j->second;
				 }

				 
			 }
			 vector<pair<int, int>> swapsToTake = stitchingResult.second;
			 vector < set<pair<int, int>>> parallelizedSwaps = parallelizeSwapPathway(swapsToTake);
			 //now we add all the swaps to the result circuit
				//all these swaps are printed/written to the output
			 for (auto z = parallelizedSwaps.begin();z != parallelizedSwaps.end();z++) {
				 for (auto y = (*z).begin();y != (*z).end();y++) {
					 cout << "swp q[" << newLastMapping[y->first] << "], q[" << newLastMapping[y->second] << "]; //swapping logical "<< y->first<< " with logical "<< y->second<<"\n";
					 //performing the swap to transform our last mapping
					 int tmp = newLastMapping[y->first];
					 newLastMapping[y->first] = newLastMapping[y->second];
					 newLastMapping[y->second] = tmp;
				 }
				 
			 }

			 //adjusting the last mapping now
			 lastMapping = make_pair(newLastMapping, lastMapping.second);

			 
			 
		 }
		 //now we have a "best mapping" and we printed out the swaps needed, if any
			//we will update the result circuit
		 //updating the result circuit:
		 //we go through every gate starting at last gates and ending at any gate in the "remaining" portion of the result
			//all these gates are printed/written to the output
		 //we are printing the output for now
		 set<GateNode*> notMapped = get<2>(bestTuple);
		 //printf("CONFIRMING NOT MAPPED GATE SET!\n");
		 for (auto z = notMapped.begin();z != notMapped.end();z++) {
			 //cout << "COULD NOT MAP GATE " << (*z)->name << "with control " << (*z)->control << " and target " << (*z)->target <<"\n";
		 }
		 queue<GateNode*> printQueue;
		 set<GateNode*> closedSet;
		 for (auto z = nextGates.begin();z != nextGates.end();z++) {
			 //adding gates to a printQueue
			 printQueue.push(*z);
			 //cout << "initial add gate " << (*z)->name << " with control " << (*z)->control << " and target " << (*z)->target << "\n";
			//printing this gate and adding its children to the printQueue, only if both parents are in the closedSet 
			 
		 }



		 //building printing set, and printing as we go
		 
		 while (printQueue.size() > 0) {
			 GateNode* nextGate = printQueue.front();
			 printQueue.pop();
			 if (closedSet.find(nextGate) == closedSet.end() && notMapped.find(nextGate)==notMapped.end()) {
				 //printing the gate
				 int control = (nextGate)->control;
				 int target = (nextGate)->target;
				 //cout << "Printing valid gate with name " << nextGate->name << " with control " << nextGate->control << "and target " << nextGate->target;
				 if (control == -1) {
					//printing out a single qubit gate
					 cout << (nextGate)->name << " q[" << get<0>(bestTuple)[target] << "]; // "<<" formerly "<< (nextGate)->name << " q["<< target<<"]\n";
				 }
				 else {
					//printing out a two qubit gate
					 cout << (nextGate)->name << " q[" << get<0>(bestTuple)[control] << "], q[" << get<0>(bestTuple)[target] << "]; //" << "formerly "<< (nextGate)->name << " q["<< control <<"], q["<< target<<"] \n";
				 }
				 //adding gate to closed set
				 closedSet.insert(nextGate);
				 //adding children to printQueue, only if both parents are in closedSet
				 if ((nextGate)->controlChild != NULL) {
					 if (nextGate->controlChild->targetParent == NULL || closedSet.find((nextGate)->controlChild->targetParent) != closedSet.end() || (closedSet.find(nextGate->controlChild->targetParent)==closedSet.end() && nextGates.find(nextGate->controlChild->targetParent)==nextGates.end())) {
						 //we add the gate to the queue
						 printQueue.push((nextGate)->controlChild);
					}
				 }
				 
				 
				 
				 if ((nextGate)->targetChild != NULL) {
					 if (nextGate->targetChild->controlParent == NULL || closedSet.find((nextGate)->targetChild->controlParent) != closedSet.end() || (closedSet.find(nextGate->targetChild->controlParent)==closedSet.end() && nextGates.find(nextGate->targetChild->controlParent)==nextGates.end())) {
						 //we add the gate to the queue
						 printQueue.push((nextGate)->targetChild);
					}
				 }
				 
			}
		 }
		 //checking if we need to account for swaps
		 nextGates = get<2>(bestTuple);
		 count++;
	 }

	 //printing our initial mapping
	 cout << "INITIAL MAPPING:\n";
	 for (auto i = initialMapping.begin();i != initialMapping.end();i++) {
		 cout << "Logical qubit q" << i->first << " mapped to physical qubit Q" << i->second << " \n";
	 }
 }

//swap optimal version of the lazy circuit builder
	//this is swap optimal
 void lazySwapOptimalCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates, vector<vector<int>> distanceMatrix) {
	 
	
	//1) get first and second layer of maximal mappings
    //2) stitch the two mappings with the lowest cost (curr Swaps in pathway + estimated swaps to stitch)
	// 3) keep repeating until we hit the last layer of maximal mappings
		//our solution here is guaranteed to be optimal


	priority_queue<SecondPriorityData> layerMappingQueue;

	//getting first layer and second layer first
		//could be the case that mapping is perfect to start, we need to check for this

	//after getting first layer and second layer, we can start pushing possibilities to the priority queue and start that logic
	vector<tuple<map<int,int>,set<int>, set<GateNode*>>> firstLayer =maximalMapper(architectureEdges, startGates);

	if (get<2>(firstLayer[0]).size() == 0) {
		//THEN WE FOUND A PERFECT MAPPING!

	}
	
	for (auto i = firstLayer.begin();i != firstLayer.end();i++) {
		vector<set<GateNode*>> endings;
		endings.push_back(get<2>(*i));
		//calling for the second layer and enqueing possibilities with cost to the priority queue
		vector<tuple<map<int, int>, set<int>, set<GateNode*>>> secondLayer = maximalMapper(architectureEdges, get<2>(*i));
		//getting cost and enqueuing
		for (auto j = secondLayer.begin();j != secondLayer.end();j++) {
			
			int cost = partialMappingStitchingCost(get<0>(*i), get<0>(*j), distanceMatrix, architectureEdges);
			vector < vector<set<pair<int, int>>>> emptySwaps;
			layerMappingQueue.push(SecondPriorityData(cost, emptySwaps,get<0>(*i),get<0>(*j),get<0>(*i), get<2>(*j),endings));
		}
	}


	
	while (layerMappingQueue.size() > 0) {
		SecondPriorityData toPursue = layerMappingQueue.top();
		layerMappingQueue.pop();

		// this is the lowest cost thing so far

		//process: 
				//1: stitch mappingOne and mappingTwo together
				//2: update swap List
				//3: update transformed mapping and initial mapping
				//4a: if end gates is empty, then we hit the end of the circuit and we can print out this swap pathway
				//4b: otherwise get the next maximal partition and enqueue the nodes with their cost
	}

	//we found the optimal swap solution: logic here is for printing it
}

 //smarter function to insert swaps (considers delay and latency) 
 void optimalLatencyOrientedCircuitBuilder() {

 }



 //function to test my maximal mapper
	//it should return mappings for every maximal portion of the circuit that can be executed without swaps
 queue<map<int, int>> getMaximalMappings(set<GateNode*> startGates,map<int,set<int>> architectureEdges) {
	 set<GateNode*> nextGates = startGates;

	 queue<map<int, int>> resultQueue;
	 int count = 1;
	 do {
		 ////printf("ON ITERATION %d\n", count);
		vector <tuple<map<int, int>, set<int>, set<GateNode*>> > result = maximalMapper(architectureEdges,nextGates);
		////printf("FINISHED MAPPING THIS SEGMENT WITH NUM MAPPINGS %u!\n",result.size());
		map<int, int> mapping = get<0>(result[0]);
		nextGates = get<2>(result[0]);
		resultQueue.push(mapping);
		count++;
	 } while (nextGates.size() > 0);

	 //now resultQueue holds all the maximal mappings for the entire circuit, in order
	 //cehcking to see if we satisfied all gates
	 return resultQueue;
 }

 //just a function to return whether a transformed circuit is a valid circuit to run
	//returns a tuple with total # of swaps, # of cycles needed to run, which are both negative if the circuit is invalid 
 pair<int,int> isValidCircuit(map<int,int> initialMapping,set<GateNode*> resultCircuit,map<int,set<int>> architectureEdges) {
	 set<GateNode*> currentLevel = resultCircuit;
	 map<int, int> currentMapping = initialMapping;
	 int numCycles=0;
	 int numSwaps = 0;
	 while (currentLevel.size() > 0) {
		//for any swaps in a level of the dependence graph, we handle those after the non-swap gates are completed 	
		 while (true) {
			 set<GateNode*>swaps;
			 for (auto i = currentLevel.begin();i != currentLevel.end();i++) {
				 GateNode* toConsider = *i;
				 string swap = "swp";
				 //we have to check if the logical qubits to swap are next to each other
					//otherwise this is not possible!
				 if (!strcmp((toConsider->name).c_str(), swap.c_str())) {
					 //then we push this gate to our swaps set
					 swaps.insert(toConsider);
				 }
				 else {
					 //check if the gate is satisfied with the current mapping
					 if (toConsider->control != -1) {
						 //we do not consider 1 qubit gates
						 if (architectureEdges[currentMapping[toConsider->control]].find(currentMapping[toConsider->target]) == architectureEdges[currentMapping[toConsider->control]].end()) {
							//then the circuit is not satisfied!	
							 return make_pair(-1,-1);
						 }
					 }
				 }
			 }
			 numSwaps += swaps.size();
			 //now going through all the swaps
			 for (auto j = swaps.begin();j != swaps.end();j++) {
				//transforming the mapping based on the swap (also seeing if the swap is legal or not)
				 GateNode* toConsider = *j;
				 if (architectureEdges[currentMapping[toConsider->control]].find(currentMapping[toConsider->target]) == architectureEdges[currentMapping[toConsider->control]].end()) {
							//then the swap is not possible!	
							 return make_pair(-1,-1);
				  }

				 //swap is possible so now we transform our initial mapping
				 int temp = initialMapping[toConsider->control];
				 initialMapping[toConsider->control] = initialMapping[toConsider->target];
				 initialMapping[toConsider->target] = temp;
		     }
			 //resetting currentLevel
				//getting the next level of the dependence graph
			 currentLevel = getNextLevel(currentLevel);
			 numCycles++;
		 }
	 }

	 return make_pair(numCycles,numSwaps);
}

 
