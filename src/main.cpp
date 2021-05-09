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
 int partialMappingStitchingCost(map<int, int> firstMap, map<int, int> secondMap, map<pair<int,int>,int> distances, map<int, set<int>> architectureEdges);
 int partialMappingCost(map<int, int> mapping, set<GateNode*> remainingGates, map<pair<int,int>,int> distances, map<int, set<int>> architectureEdges);
 pair<map<int, int>, vector<pair<int, int>>> stitchMappings(pair<map<int, int>, set<int>> mapOne, pair<map<int, int>, set<int>> mapTwo, map<int, set<int>> architectureEdges,map<pair<int,int>,int> distances);
//void lazySwapCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates, vector<vector<int>> distanceMatrix);
void betterLazySwapCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates);
void optimalLazySwapCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates);
map<int, set<int>> getCroppedArchitecture(map<int, int> filledMapping, map<int, set<int>> architectureEdges);
vector<map<int, int>> fillPartialMapping(map<int, int> mapping, map<int, set<int>> architectureEdges);
map<int, int> pickBestFilledMapping(map<int, int> mapping, map<int, set<int>> architectureEdges, set<GateNode*> remainingGates, map<pair<int,int>,int> distances);
map<pair<int, int>, int> getDistances(map<int, set<int>> architectureEdges);
int getNaiveDistance(map<int, set<int>> architectureEdges, int Q1, int Q2);
vector<vector<pair<int, int>>> getSwapPathways(map<int, int> one, map<int, int> two, map<pair<int, int>, int> distances, map<int, set<int>> architectureEdges);


//global distance matrix
 //vector<vector<int>> distanceMatrix;

 struct PriorityData {
	 int cost;
	 vector<pair<int, int>> swapsTaken;
	 //map<int, int> initialMapping;
	 map<int, int> transformedMapping;

	 PriorityData(int cost, vector<pair<int, int>> swapsTaken, map<int, int> transformedMapping) {
		 this->cost = cost;
		 this->swapsTaken = swapsTaken;
		 //this.initialMapping = third;
		 this->transformedMapping = transformedMapping;
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
	 vector<vector<pair<int,int>>> swapsTaken;
	 //keeping track of the initial mapping for the entire process
	 map<int, int> initialMapping;
	 //this is the next proposed mapping in the next layer of maximal mappings to stitch to
	 map<int, int>  targetMapping;
	//this is our current transformedMapping in whatever maximal mapping layer we are at
	 map<int, int> transformedMapping;
	//keeping track of which gate set we are at 
	 set<GateNode*> firstGates;
	 //keeping track of each section of end gates for final printing
	 vector<queue<GateNode*>> endingSections;

	 SecondPriorityData(int heuristicCost,  vector<vector<pair<int,int>>> maximalPartitionSwaps, map<int, int> initialMap, map<int, int> targetMap, map<int,int> transformedMap, set<GateNode*> currentGates, vector<queue<GateNode*>> endings) {
		 cost = heuristicCost;
		 swapsTaken = maximalPartitionSwaps;
		 initialMapping = initialMap;
		 targetMapping = targetMap;
		 transformedMapping = transformedMap;
		 firstGates = currentGates;
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
	//distanceMatrix = getPhysicalDistancesFloydWarshall(architectureEdges);

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
	//betterLazySwapCircuitBuilder(architectureEdges, firstGates);
	/*optimal test*/
	optimalLazySwapCircuitBuilder(architectureEdges, firstGates);
	//Exit the program:
	return 0;
}

//maps each pair of qubits to a distance value
	map<pair<int, int>, int> getDistances(map<int,set<int>> architectureEdges) {
		map<pair<int, int>, int> distanceMap;
		for (auto i = architectureEdges.begin();i != architectureEdges.end();i++) {
			for (auto j = architectureEdges.begin();j != architectureEdges.end();j++) {
				int distance = getNaiveDistance(architectureEdges, i->first, j->first);
				//updating distance value

				distanceMap[make_pair(i->first, j->first)] = distance;
			}
		}

		//printf("PRINTING DISTANCES!\n");
		//for (auto i = distanceMap.begin();i != distanceMap.end();i++) {
		//	printf("%d to %d has distance %d\n", i->first.first, i->first.second, i->second);
		//}
		return distanceMap;
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
	//if fast is true, the algorithm will try and rule out partial mappings at every step of the propogation process
		//this may result in a maximal mapping that is not actually maximal, but it is faster for larger circuits
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
	//set<GateNode*> otherClosedSet; 
	
	//we do the first iteration
	GateNode* firstGate = mappingQueue.front();	
	mappingQueue.pop_front();
	//otherClosedSet.insert(firstGate);
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
	printf("GATES TO MAP SIZE: %lu\n", mappingQueue.size());
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
		
		//maybe do fast optimization here
		//we should cut down on mappings to propogate in order to save on execution time
		//this below step is just for the large benchmarks, without this, they balloon and use too much memory 
		while (maximalMappingQueue.size() > 60000) {
			//throwing out any element
			maximalMappingQueue.pop();
		}

		unsigned int internalStep =0;
		unsigned int size = maximalMappingQueue.size();
		//cout << "CURRENT MAXIMAL MAPPING QUEUE SIZE: " << size << "\n";
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
		//otherClosedSet.insert(toSatisfy);
		//updating counter of number of mappings in best mapping so far
		lastMapped = internalStep;
	}
	//printf("PRINTING GATES CONSIDERED IN THIS ROUND OF MAPPING:\n");
	//for (auto z = otherClosedSet.begin();z != otherClosedSet.end();z++) {
		//cout << "Gate name: " + (*z)->name << " control: " << (*z)->control << " target: " << (*z)->target <<"\n";
	//}
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
		/*
		for (auto j = otherClosedSet.begin();j != otherClosedSet.end();j++) {
			if (completedGates.find((*j)) == completedGates.end()) {
				//then this wasnt mapped
				//we need to see if its parents were mapped or not
				finalSet.insert(*j);
			}
		}
		*/
		//i just put the couldntMap gates and children of completedGates into the nextGates set
		for (auto j = completedGates.begin();j != completedGates.end();j++) {
			if ((*j)->controlChild != NULL) {
				//checking if child is not in completed gates or couldnt map set
				if (completedGates.find((*j)->controlChild) == completedGates.end() && couldntMap.find((*j)->controlChild) == couldntMap.end()) {
					finalSet.insert((*j)->controlChild);
				}
			}

			if ((*j)->targetChild != NULL && (*j)->targetChild!=(*j)->controlChild) {
				if (completedGates.find((*j)->targetChild) == completedGates.end() && couldntMap.find((*j)->targetChild) == couldntMap.end()) {
					finalSet.insert((*j)->targetChild);
				}
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
	//cout << "Gates Considered Size: " << otherClosedSet.size() << "\n";
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
	 printf("STARTING SETUP!\n");
	 int numEdges = 0;
	 for (auto i = architectureEdges.begin();i != architectureEdges.end();i++) {
		 numEdges += (i->second).size();
	 }
	 vector<vector<int>> distances(architectureEdges.size(), vector<int>(architectureEdges.size(),numEdges));
	 printf("FINISHED SETUP!\n");

	 //setting up the algorithm
	 for (auto i = architectureEdges.begin();i != architectureEdges.end();i++) {
		 distances[i->first][i->first] = 0;
		 for (auto j = (i->second).begin();j != (i->second).end();j++) {
			 //distance 1 edges
			 distances[i->first][*j] = 1;
		 }
	 }
	 printf("FINISHED MORE SETUP!\n");
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
 int partialMappingCost(map<int,int> mapping, set<GateNode*> remainingGates, map<pair<int,int>,int> distances, map<int,set<int>> architectureEdges) {
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
			 cost += distances.find(make_pair(mapping[toAppraise->control], mapping[toAppraise->target]))->second;
			 //cost += getDistance(distanceMatrix, architectureEdges ,mapping[toAppraise->control], mapping[toAppraise->target]) - 1;
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
 int partialMappingStitchingCost(map<int,int> firstMap, map<int,int> secondMap,map<pair<int,int>,int> distances, map<int,set<int>> architectureEdges) {
	//all we do is get the maximum distance between any 2 mapped logical qubits in both mappings	
	 int maxDistance = 0;
	 auto iter = firstMap.begin();
	 while (iter != firstMap.end()) {
		 if (secondMap.find(iter->first) != secondMap.end()) {
			 //then this logical qubit (iter->first) is in both partial mappings
			 int calculatedDistance = distances.find(make_pair(iter->second, iter->first))->second;
			 //int calculatedDistance = getDistance(distanceMatrix, architectureEdges, firstMap[iter->first], secondMap[iter->first]);
			 if ( calculatedDistance> maxDistance) {
				 maxDistance = calculatedDistance;
			 }
		}
		 iter++;
	 }
	 //now we have the minimum # (lower bound) on swaps needed in order to transform a mapping into another mapping
	 return maxDistance;

 }

 //will enumerate possibilities for remainder of partial mapping
 vector<map<int, int>> fillPartialMapping(map<int,int> mapping, map<int,set<int>> architectureEdges) {
	 vector<map<int, int>> results;
	 queue<map<int, int>> mapFillQueue;
	 

	 //initially enqueueing map and unmapped set
	 mapFillQueue.push(mapping);
	 while (mapFillQueue.size() > 0) {
		 map<int, int> state = mapFillQueue.front();
		 mapFillQueue.pop();
		 
		 if (state.size() == numLogical) {
			 //then we mapped all the unmapped qubits in this state
			 results.push_back(state);
			 continue;
		 }

		 //finding an unmapped neighbor
		 for (auto j = state.begin();j != state.end();j++) {
			 //keeping track of stuff we considered, so we do not consider it again
			 set<int> considered;
			 for (auto z = architectureEdges[j->second].begin();z != architectureEdges[j->second].end();z++) {
				 bool found = false;
				 for (auto i = state.begin();i != state.end();i++) {
					 if (i->second == *z) {
						 found = true;
						 break;
					 }
				 }
				 if (!found && considered.find(*z)==considered.end()) {
					 //then this is a possible mapping to fill
					 considered.insert(*z);
					 for (int k = 0;k <numLogical;k++) {
						 if (state.find(k) == state.end()) {
							map<int, int> newState = state;
							newState[k] = *z;
							//pushing this state back to queue
							mapFillQueue.push(newState);
						 }
						 
					 }
				 }
			 }
		 }
	 }
	 return results;
 }

 //for the fast version of my algorithm
	//goes through all the mappings, and picks 
 map<int, int> pickBestFilledMapping(map<int,int> mapping, map<int,set<int>> architectureEdges,set<GateNode*> remainingGates,map<pair<int,int>,int> distances) {
	 vector<map<int, int>> filledMappings = fillPartialMapping(mapping, architectureEdges);
	 map<int, int> bestMap;
	 int bestMapCost = 0;
	 for (auto i = filledMappings.begin();i != filledMappings.end();i++) {
		 int cost = partialMappingCost(*i, remainingGates, distances, architectureEdges);
		 if (bestMap.size() == 0) {
			 bestMap = *i;
			 bestMapCost = cost;
		 }
		 else if (cost < bestMapCost) {
			 bestMap = *i;
			 bestMapCost = cost;
		 }
	 }
	 return bestMap;
 }


 //optimal solver

 void optimalSwapSolver(map<int, set<int>> architectureEdges, set<GateNode*> startGates) {
	 

 }

 

 

 

 //function to get a list of swaps needed to transform one mapping to another
	//returns as the first element of the pair the modified initial mapping (additional logical qubits mapped to mapOne, if necessary)
	//second element of the returned pair is the list of swaps to take
	//if fast is true, then optimization on cutting down priority_queue size is taken
 pair<map<int, int>, vector<pair<int, int>>> stitchMappings(pair<map<int, int>, set<int>> mapOne, pair<map<int, int>, set<int>> mapTwo, map<int, set<int>> architectureEdges, map<pair<int, int>, int> distances) {
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

	 map<int, int> firstMap = mapOne.first;
	 map<int, int> transformedMap = firstMap;
	 set<int> firstArchMap = mapOne.second;
	 map<int, int> secondMap = mapTwo.first;
	 set<int> secondArchMap = mapTwo.second;


	 //enqueueing the base state
	 printf("ENQUEING BASE STATE!\n");
	 swapQueue.push(PriorityData(partialMappingStitchingCost(firstMap, secondMap, distances, architectureEdges), vector<pair<int, int>>(), transformedMap));
	 int iterCount = 0;
	 while (swapQueue.size() > 0) {
		 printf("SIZE OF SWAP QUEUE: %d\n", swapQueue.size());

		 while (swapQueue.size() > 60000 ) {
			 //just an optimization to cut down on mappings "we think" are bad, especially useful for large benchmarks
			 //we will only keep the top 3 mappings from swapQueue
			 priority_queue<PriorityData> tmp;
			 for (int i = 0;i < 3;i++) {
				 tmp.push(swapQueue.top());
				 swapQueue.pop();
			 }
			 swapQueue = tmp;
		 }
		  //printf("GOT INSIDE!\n");
		 //getting the lowest cost path from the queue
		 PriorityData result = swapQueue.top();
		 //tuple<int,vector<pair<int, int>>,map<int, int>,map<int, int>>  result = swapQueue.top();
		 swapQueue.pop();
		 //vector<pair<int, int>> swapsTaken = get<1>(result);
		 //map<int, int> initialMapping = get<2>(result);
		// map<int, int> transformedMapping = get<2>(result);
		 //int currCost = get<0>(result);
		 vector<pair<int, int>> swapsTaken = result.swapsTaken;
		 //map<int, int> initialMapping = result.initialMapping;
		 map<int, int> transformedMapping = result.transformedMapping;
		 
		 int currCost = result.cost;
		 //printf("CURRENT TRANSFORMED MAPPING:\n");
		 //for (auto o = transformedMapping.begin();o != transformedMapping.end();o++) {
			 //printf("%d mapped to %d, ", o->first, o->second);
		 //}
		 if (swapsTaken.size() != 0) {
			 //printf("LAST CHOSEN SWAP WAS BETWEEN %d and %d with Cost %d\n", swapsTaken[swapsTaken.size() - 1].first, swapsTaken[swapsTaken.size() - 1].second, currCost);
		 }
		 //printf("WITH COST: %d", currCost);
		 //printf("\n");
		 //debug

		 //checking equality 
		 bool equal = true;
		 for (auto h = secondMap.begin();h != secondMap.end();h++) {

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
			 return make_pair(firstMap, swapsTaken);
		 }


		 //else we need to think in terms of 1 swap granularities
		 //we have some cases to think about:
			//1) suppose a qubit exists in mapping 1 but not in mapping 2
				//this is no issue because then this qubit does not contribute to the next section of the dependence graph
			//2) suppose a logical qubit is mapped in mapping 2 but not in mapping 1
				//then that logical qubit must be involved in an adjacent cnot gate
				//therefore, we try mapping common qubits first, and as we consider swaps with unmapped neighbors, we assign qubits
			//3) if a logical qubit is mapped in 1 and 2
				//then we just consider possible swaps of distance 1 away that are close
		 for (auto j = transformedMapping.begin();j != transformedMapping.end();j++) {
			 if (secondMap.find(j->first) == secondMap.end()) {
				 //then this qubit does not contribute to the next maximal mapping 
				 //printf("LOGICAL QUBIT %d DOES NOT CONTRIBUTE TO NEXT MAPPING!\n", j->first);

				 continue;
			 }
			 if (transformedMapping[j->first] != secondMap[j->first]) {
				 pair<int, int> bestSwap = make_pair(-1, -1);
				 int distance = (distances.find(make_pair(j->second, secondMap[j->first])))->second;
				 for (auto i = transformedMapping.begin();i != transformedMapping.end();i++) {
					 if (architectureEdges[j->second].find(i->second) != architectureEdges[j->second].end()) {
						 //then this is a possible swap
						 int firstDistance = (distances.find(make_pair(i->second, secondMap[j->first])))->second;
						 int secondDistance = (distances.find(make_pair(j->second, secondMap[j->first])))->second;

						 if (firstDistance < distance) {
							 //updating swap and distance
							 bestSwap = make_pair(i->first, j->first);
							 distance = firstDistance;
						 }
					 }
				 }

				 //doing the best Swap
				 vector<pair<int, int>> newSwaps = swapsTaken;
				 newSwaps.push_back(bestSwap);
				 //printf("CONSIDERING SWAPPING %d and %d\n", i->first, j->first);

				 //updating transformed mapping
				 map<int, int> updatedMapping = transformedMapping;
				 int tmp = updatedMapping[bestSwap.first];
				 updatedMapping[bestSwap.first] = updatedMapping[bestSwap.second];
				 updatedMapping[bestSwap.second] = tmp;

				 //updating cost
				 int cost = newSwaps.size() + partialMappingStitchingCost(updatedMapping, secondMap, distances, architectureEdges);

				 //enqueueing
				 swapQueue.push(PriorityData(cost, newSwaps, updatedMapping));
			 }
		 }
	 }
	 return make_pair(map<int, int>(), vector<pair<int, int>>());
 }
					/*
						 
						 if (firstDistance < secondDistance){
							 //possible swap

							 //checking redundant
							 if (swapsTaken.size() > 0) {
								 pair<int, int> lastSwap = swapsTaken[swapsTaken.size() - 1];
								 if ((lastSwap.first == i->first && lastSwap.second == j->first) || (lastSwap.first == j->first && lastSwap.second == i->first)) {
									 //redundant swap
									 continue;
								 }
							 }
							 //not a redundant swap and we can enqueue
							 vector<pair<int, int>> newSwaps = swapsTaken;
							 newSwaps.push_back(make_pair(i->first, j->first));
							 //printf("CONSIDERING SWAPPING %d and %d\n", i->first, j->first);

							 //updating transformed mapping
							 map<int, int> updatedMapping = transformedMapping;
							 int tmp = updatedMapping[i->first];
							 updatedMapping[i->first] = updatedMapping[j->first];
							 updatedMapping[j->first] = tmp;

							 //updating cost
							 int cost = newSwaps.size() + partialMappingStitchingCost(updatedMapping, secondMap,distances, architectureEdges);

							 //enqueueing
							 swapQueue.push(PriorityData(cost, newSwaps, updatedMapping));
							 

							 //for fast mapper, we will only consider 1 swap at a time in the queue, or it baloons
							 break;

						 }
					 }
				 }
			 }
		 }
	 }
	 return make_pair(map<int, int>(), vector<pair<int, int>>());
 }
 */
				 /*
				 //printf("CASE 1 NEED TO SWAP!\n");
				 //will need to consider 1 swap applied to j->first in order to move it closer to secondMap[j->first]
					//of course this 1 swap needs to bring j->first closer to secondMap[j->first] or vice versa
				 // for loop below considering a single swap applied to logical qubit mapping j->first
				 for (auto i = architectureEdges[transformedMapping[j->first]].begin();i != architectureEdges[transformedMapping[j->first]].end();i++) {
					 //need to check if this qubit is mapped to in the current mapping or not
						//if so, we consider a swap with the associated logical qubit
						//if not, then we check if there are qubits that are mapped in the result that are not mapped in current
							//this is a possible assignment for them
					 //is this physical qubit closer to secondMap[j->first] or not? 
						//we only consider the swap if this is the case
					 /*
					 if (distanceMatrix[*i][secondMap[j->first]] >= distanceMatrix[transformedMapping[j->first]][secondMap[j->first]]) {
						 //swap is useless
						 continue;
					 }
					 */
		/*
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
						 swapQueue.push(PriorityData(newCost, newSwaps, newTransformedMapping));
					 	  
			 }
		 }
		 iterCount++;
	 }
	 return make_pair(map<int,int>(), vector<pair<int,int>>());
 }
*/
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

 //gets a cropped architecture from a filled mapping
	//this will help for future 
 map<int, set<int>> getCroppedArchitecture(map<int,int> filledMapping, map<int, set<int>> architectureEdges) {
	 map<int, set<int>> newArchitecture;
	 for (auto i = filledMapping.begin();i != filledMapping.end();i++) {
		 //adding the mapped physical qubit to the architecture, and all its mapped neighbors
			//since this is a full mapping, at least 1 exists
		 newArchitecture[i->second] = set<int>();
		 for (auto z = filledMapping.begin();z != filledMapping.end();z++) {
			 if (architectureEdges[i->second].find(z->second) != architectureEdges[i->second].end()) {
				 //then we add this edge
				 newArchitecture[i->second].insert(z->second);
			 }
		 }
		 
		 if (newArchitecture[i->second].size() == 0) {
			 printf("SOMETHING WENT WRONG HERE!\n");
		 }
	 }
	 printf("PRINTING CROPPED ARCHITECTURE!\n");
	 for (auto z = newArchitecture.begin();z != newArchitecture.end();z++) {
		 printf("PHYSICAL QUBIT %d mappings: ", z->first);
		 for (auto y = (z->second).begin();y != (z->second).end();y++) {
			 printf(" %d, ", *y);
		 }
		 printf("\n");
	 }
	 //returning this new architecture
	 return newArchitecture;
 }

 //helper function to actually schedule the transformed circuit given maximally mapped gates, swaps, and an initial mapping and print the results
	//also prints the initial mapping,cycle depth, and number of gates as requested
	//this takes latencies into account and schedules the gates according to the latency
	//interweave is true if swaps should try and be interleaved in a maximal partition
	//if interweave is false, swaps are just tacked on to the end of a maximal partition
 void scheduleCircuit(vector<queue<GateNode*>> gateMappings,vector<vector<pair<int,int>>> levelSwaps,map<int,int> initialMapping, int latencySingle, int latencyDouble, int latencySwap,map<int,set<int>> architectureEdges) {
	//idea, we just schedule level by level, and whatever is not finished after every cycle, we push to the next	
		//this is the trivial scheduling for the non-interleaved case
	 int numCycles = 0;
	 int numGates = 0;
	 GateNode* dummy = new GateNode();
	 map<int, int> toTransform = initialMapping;
	 vector<pair<GateNode*, int>> logicalSchedule; 
	 for (int i = 0;i < numLogical;i++) {
		 //pushing initial state of qubits in circuit
		 logicalSchedule.push_back(make_pair(dummy, -1));
	 }
	 //implementation of greedy scheduling  
	 for (int i =0;i<gateMappings.size();i++) {
		//putting queue into a set
		 queue<GateNode*> satisfied = gateMappings[i];
		 set<GateNode*> toSchedule;
		 set<GateNode*> inProgress;
		 set<GateNode*> finished;
		 while (satisfied.size() > 0) {
			 GateNode* toAdd = satisfied.front();
			 satisfied.pop();
			 toSchedule.insert(toAdd);
		 }

		 vector<set<pair<int, int>>> parallelizedSwaps;

		 if (i != gateMappings.size() - 1) {
			 parallelizedSwaps = parallelizeSwapPathway(levelSwaps[i]);
		 }

		 //debug
		 printf("PRINTING ALL GATES TO SCHEDULE!\n");
		 for (auto j = toSchedule.begin();j != toSchedule.end();j++) {
			 cout << (*j)->name << " between control: " << (*j)->control << " and target: " << (*j)->target << " \n";
		 }
				 
		 //starting the scheduling
		 int numTotalGates = toSchedule.size();
		 int numSwapsAdded = 0;
		 while (finished.size() != numTotalGates + numSwapsAdded) {
			 for (auto j = toSchedule.begin();j != toSchedule.end();) {
				 GateNode* toConsider = *j;
				 if (toConsider->control == -1) {
					 //single qubit gate to schedule

					  //checking if parents are done
					 if (toConsider->targetParent == NULL ||
						 finished.find(toConsider->targetParent) != finished.end() ||
						 (toSchedule.find(toConsider->targetParent) == toSchedule.end() && inProgress.find(toConsider->targetParent) == inProgress.end())) {
						 //then the parent is either done or not considered, so we are good with scheduling
						 //need to check if this qubit is busy or not
						 if (logicalSchedule[toConsider->target].second == -1) {
							 //we can schedule
							 logicalSchedule[toConsider->target] = make_pair(toConsider, latencySingle);
							 inProgress.insert(toConsider);
							 //removing gate from toSchedule
							 j = toSchedule.erase(j);

							 //incrementing the gate number
							 numGates++;
						 }
						 else {
							 //moving on
							 j++;
						 }
					 }
					 else {
						 //moving on
						 j++;
					 }
				 }
				 else {
					 //two qubit gate
					 //need to check if both parents are satisfied
					 if ((toConsider->targetParent == NULL ||
						 finished.find(toConsider->targetParent) != finished.end() ||
						 (toSchedule.find(toConsider->targetParent) == toSchedule.end() && inProgress.find(toConsider->targetParent) == inProgress.end())) &&
						 ((toConsider->controlParent == NULL) || finished.find(toConsider->controlParent) != finished.end() ||
							 (toSchedule.find(toConsider->targetParent) == toSchedule.end() && inProgress.find(toConsider->targetParent) == inProgress.end()))) {
						 //then both parents are satisfied, lets try to map these qubits
						 if (logicalSchedule[toConsider->target].second == -1 && logicalSchedule[toConsider->control].second == -1) {
							 //then we can schedule this gate for both qubits
							 logicalSchedule[toConsider->target] = make_pair(toConsider, latencyDouble);
							//removing this gate from schedule list
							 inProgress.insert(toConsider);
							 //removing gate from toSchedule
							 j = toSchedule.erase(j);
							 //incrementing gate number
							 numGates++;

							 //debug--> making sure that this is valid
							 bool found = false;
							 int firstPhys = toTransform[toConsider->control];
							 int secondPhys = toTransform[toConsider->target];
							 for (auto z = architectureEdges[firstPhys].begin();z != architectureEdges[firstPhys].end();z++) {
								 if (*z == secondPhys) {
									 found = true;
									 break;
								 }
							 }
							 if (!found) {
								 cout << " SOMETHING WENT TERRIBLY WRONG IN SCHEDULING A QUBIT!\n";
							 }
						 }
						 else {
							 //moving on 
							 j++;
						 }
					 }
					 else {
						 //moving on 
						 j++;
					 }
				 }
			 }

			 //we check here if we can put in a swap or not
				//if a space is empty and the two qubits are not in use and they will not be in use for the rest, then we can safely insert the swap
			 printf("STARTED SWAP CONSIDERATION!\n");
			 if (i != gateMappings.size() - 1) {
				
				 //then we can consider swaps
					 //we already have swaps separated by parallel swaps
					 //so we try and do all the swaps in a section first, and then move on to the next set of swaps and so on
					 //special case is if the swaps in the current set do not conflict with the swaps in the next set at any time
						 //if so, we can mash them together!			
						 //specifically, after swapping, there might be some swaps in the second set that do not conflict with the current set
							 //if so, we can add them, and adjust the entire set while we are at it
				 for (auto z = parallelizedSwaps.begin();z != parallelizedSwaps.end();) {
					 for (auto y = (*z).begin();y != (*z).end();) {
						 //seeing if we can add a swap or not
						 if (logicalSchedule[y->first].second == -1 && logicalSchedule[y->second].second == -1) {
							 //then we go through all the next gates to map and see if these 2 qubits are involved: if not, then we can go ahead	
							 bool anyQubitInvolved = false;
							 for (auto p = toSchedule.begin();p != toSchedule.end();p++) {
								 if ((*p)->target == y->first || (*p)->target == y->second || (*p)->control == y->first || (*p)->target == y->second) {
									 //then we shouldnt swap
									 anyQubitInvolved = true;
									 break;
								 }
							 }
							 if (anyQubitInvolved) {
								 //then we cant swap this 
								 y++;
							 }
							 else {
								 //then we can swap these qubits
								 printf("CULPRIT!\n");
								 GateNode* newSwapGate = new(GateNode);
								 newSwapGate->name = "swp";
								 newSwapGate->target = y->second;
								 newSwapGate->control = y->first;
								 newSwapGate->targetParent = NULL;
								 newSwapGate->controlParent = NULL;
								 printf("AHA!\n");

								 //adding this gate to the schedule and removing it from this set
								 logicalSchedule[y->first] = make_pair(newSwapGate, latencySwap);
								 y = (*z).erase(y);

								 //incrementing number of gates
								 numGates++;
								 //incrementing number of swaps added
								 numSwapsAdded++;

								 printf("MAYBE NOT?!\n");
							 }
						 }
						 else {
							 //then we cant add the swap
							 y++;
						 }
					 }
					 //checking if the size is zero, for now, if size is not zero, we break
					 //if size is zero, we can remove this set and continue to the next set of swaps
					 if ((*z).size() > 0) {
						 //cant do anything else
						 break;
					 }
					 else {
						 //we can remove this and move on
						 printf("DID A REMOVE!\n");
						 z = parallelizedSwaps.erase(z);
					 }
				 }
			 }

			 printf("GOT HERE!\n");

			 //we have scheduled all that we can right now, let us move onto the decrement step for cycle logic
				//we decrement all in progress gates by 1
			 for (auto z = logicalSchedule.begin();z != logicalSchedule.end();z++) {
				 if (z->second > 0) {
					 //then lets try decrementing this
					 GateNode* work = z->first;
					 int cyclesLeft = z->second;
					 cyclesLeft--;
					 if (cyclesLeft == 0) {
						 //then this gate is done executing, we should remove it from the inProgress set and from logicalSchedule
						 inProgress.erase(work);
						 finished.insert(work);
						 //if this gate is a swap gate, then we transform the mapping
						 *z = make_pair(dummy, -1);
						 if (work->name == "swp") {
							 //we do the swap change
							 //debug --> checking if swap is valid or not
							 int physOne = toTransform[work->control];
							 int physTwo = toTransform[work->target];
							 bool found = false;
							 for (auto k = architectureEdges[physOne].begin();k != architectureEdges[physOne].end();k++) {
								 if (*k == physTwo) {
									 found = true;
									 break;
								 }
							 }
							 if (!found) {
								 cout << "SOMETHING WENT TERRIBLY WRONG WITH MAPPING!\n";
							 }


							 //doing the swap
							 int tmp = toTransform[work->control];
							 toTransform[work->control] = toTransform[work->target];
							 toTransform[work->target] = tmp;
						 }

						 //printing the gate
						 if (work->control == -1) {
							 cout << work->name << " q[" << toTransform[work->target] << "] ; // formerly on logical qubit q[" << work->target << "]\n";
						 }
						 else {
							 cout << work->name << " q[" << toTransform[work->control] << "], q[" << toTransform[work->target] << "] ; // between logical q[" << work->control << "], q[" << work->target << "]\n";
						 }
					 }
				 }
			 }
			 //cycle is over 
			 numCycles++;
		 }
		 
	 }
	 cout << "NUMBER OF CYCLES: " << numCycles<<"\n";
	 cout << "NUMBER OF GATES: " << numGates << "\n";

	 /*
		 //basic implementation for now
		 //just tacking on swaps to the end of a portion of the mapping

		 //idea for representing scheduling, we have a vector of size numLogical of type pair<GateNode*,int> 
		 //the gateNode is the gate that is being run , and the second part of the int is the time remaining

		 //we already have a level order for each of the gates and a set of swaps after each set of gates, so no problem here
			//if a logical qubit is not doing any work, its space will have (NULL,-1)

		 //idea, greedy scheduling
				//just schedule what we can, if there is a gap and 2 qubits are not being used in the rest, we can start the swap
				//basically, use slack space here
		 map<int, int> toTransform = initialMapping;
		 vector<pair<GateNode*, int>> logicalSchedule;
		 for (int i = 0;i < numLogical;i++) {
			 logicalSchedule.push_back(make_pair(NULL, -1));
		 }
		 int iterations = 0;
		 for (auto z = gateMappings.begin();z != gateMappings.end();z++) {
			 if (iterations != 0) {
				//need to put in swaps here
				 vector<pair<int, int>> swaps = levelSwaps[i - 1];
				 vector<set<pair<int, int>>> parallelizedSwaps = parallelizeSwapPathway(swaps);

				 //size of this vector is the number of cycles we will need
				 numCycles += parallelizedSwaps.size();

				 //transforming the mapping and printing the swap gate
				 for (auto i = parallelizedSwaps.begin();i != parallelizedSwaps.end();i++) {
					 for (auto j = (*i).begin();j != (*i).end();j++) {
						 pair<int, int> swap = *j;
						 cout << "swp q[" << toTransform[swap.first] << "], q[" << toTransform[swap.second] << "]; // between logical q[" << swap.first << "], q[" << swap.second << "]\n";
						 numGates++;
						 //transforming the mapping
						 int tmp = toTransform[swap.first];
						 toTransform[swap.first] = toTransform[swap.second];
						 toTransform[swap.second] = tmp;
					}
				 }
			 }
			 queue<GateNode*> mappings = *z;
			 while (mappings.size() > 0) {
				 //fit whatever we can to the logical qubits
				 //once we hit a gate that cannot be mapped to a logical qubit, thats the end of the level mapping
				 GateNode* toSchedule = mappings.front();
				 //see if this can be assigned to the qubit or not
				 int control = toSchedule->control;
				 int target = toSchedule->target;
				 bool couldntSchedule = false;
				 if (control == -1) {
					 //single qubit gate
					 if (logicalSchedule[target].first == NULL) {
						 //placing this in the schedule
						 //1 Latency for now in code, will change
						 pair<GateNode*, int> toInsert = make_pair(toSchedule, 1);
						 //removing the gate from mappings because it is scheduled
						 mappings.pop();
					 }
					 else {
						 couldntSchedule = true;
					 }
				 }
				 else {
					 //2 qubit gate
					 if (logicalSchedule[target].first == NULL && logicalSchedule[control].first == NULL) {
						 //placing this in the schedule
						 //1 latency in code for now, will change
						 pair<GateNode*, int> toInsert = make_pair(toSchedule, 1);
						 //removing the gate from mappings because it is scheduled
						 mappings.pop();
					 }
					 else {
						 couldntSchedule = true;
					 }

				 }

				 if (couldntSchedule) {
					 //done scheduling for now, we simulate a cycle
					 for (auto j = logicalSchedule.begin();j != logicalSchedule.end();j++) {
						 if (j->second > 0) {
							 //decrement
							 int newCyclesLeft = j->second - 1;
							 if (newCyclesLeft == 0) {
								 //then the gate is done, we can print it
								 GateNode* toPrint = j->first;
								 int control = toPrint->control;
								 int target = toPrint->target;
								 if (control == -1) {
									 //printing single qubit gate
									 cout << toPrint->name << " q[" << toTransform[target] << "]; // formerly on qubit q[" << target << "]";
									 numGates++;
									 //removing gate from schedule list
									 (*j) = make_pair(NULL, -1);
								 }
								 else {
									 //printing two qubit gate
									 cout << toPrint->name << "q[" << toTransform[control] << "], q[" << toTransform[target] << "]; // formerly between q[" << control << "], q[" << target << "] \n";
									 numGates++;
									 //removing gate from schedule list
									 logicalSchedule[control] = make_pair(NULL, -1);
									 logicalSchedule[target] = make_pair(NULL, -1);
								 }
							 }
						}
					 }
					 //we did one cycle
					 numCycles++;
				 }

			 }

			 iterations++;
		 }

	*/	
	 
 }

 vector<vector<pair<int, int>>> getSwapPathways(map<int,int> one, map<int,int> two, map<pair<int,int>,int> distances,map<int,set<int>> architectureEdges) {
	priority_queue<PriorityData> swapQueue;
	vector <vector<pair<int, int>>> results;
	 map<int, int> firstMap = one;
	 map<int, int> transformedMap = one;
	 map<int, int> secondMap = two;


	 //enqueueing the base state
	 //printf("ENQUEING BASE STATE!\n");
	 swapQueue.push(PriorityData(partialMappingStitchingCost(firstMap, secondMap, distances, architectureEdges), vector<pair<int, int>>(), transformedMap));
	 int iterCount = 0;
	 while (swapQueue.size() > 0) {
		 //printf("SIZE OF SWAP QUEUE: %d\n", swapQueue.size());

		 while (swapQueue.size() > 60000 ) {
			 //just an optimization to cut down on mappings "we think" are bad, especially useful for large benchmarks
			 //we will only keep the top 3 mappings from swapQueue
			 priority_queue<PriorityData> tmp;
			 for (int i = 0;i < 3;i++) {
				 tmp.push(swapQueue.top());
				 swapQueue.pop();
			 }
			 swapQueue = tmp;
		 }
		  //printf("GOT INSIDE!\n");
		 //getting the lowest cost path from the queue
		 PriorityData result = swapQueue.top();
		 //tuple<int,vector<pair<int, int>>,map<int, int>,map<int, int>>  result = swapQueue.top();
		 swapQueue.pop();
		 //vector<pair<int, int>> swapsTaken = get<1>(result);
		 //map<int, int> initialMapping = get<2>(result);
		// map<int, int> transformedMapping = get<2>(result);
		 //int currCost = get<0>(result);
		 vector<pair<int, int>> swapsTaken = result.swapsTaken;
		 //map<int, int> initialMapping = result.initialMapping;
		 map<int, int> transformedMapping = result.transformedMapping;
		 
		 int currCost = result.cost;
		 //printf("CURRENT TRANSFORMED MAPPING:\n");
		 //for (auto o = transformedMapping.begin();o != transformedMapping.end();o++) {
			 //printf("%d mapped to %d, ", o->first, o->second);
		 //}
		 if (swapsTaken.size() != 0) {
			 //printf("LAST CHOSEN SWAP WAS BETWEEN %d and %d with Cost %d\n", swapsTaken[swapsTaken.size() - 1].first, swapsTaken[swapsTaken.size() - 1].second, currCost);
		 }
		 //printf("WITH COST: %d", currCost);
		 //printf("\n");
		 //debug

		 //checking equality 
		 bool equal = true;
		 for (auto h = secondMap.begin();h != secondMap.end();h++) {

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
			 int bestSwaps = swapsTaken.size();
			 results.push_back(swapsTaken);
			 // WE WILL RETURN ALL SUCH BEST SWAP PATHWAYS
			 while (true) {
				 if (swapQueue.size() == 0) {
					 break;
				 }
				 PriorityData nextRes = swapQueue.top();
				 swapQueue.pop();
				 map<int, int> transformedMap = nextRes.transformedMapping;
				 bool satisfied = true;
				 for (auto z = secondMap.begin();z != secondMap.end();z++) {
					 if (z->second != transformedMap[z->first]) {
						 satisfied = false;
						 break;
					 }
				 }
				 if (satisfied && nextRes.swapsTaken.size() == bestSwaps) {
					 //push this also
					 results.push_back(nextRes.swapsTaken);
				 }
				 else {
					 //something is wrong, which means the rest of the solutions are not optimal
					 break;
				 }


			 }
			 return results;


		 }


		 //else we need to think in terms of 1 swap granularities
		 //we have some cases to think about:
			//1) suppose a qubit exists in mapping 1 but not in mapping 2
				//this is no issue because then this qubit does not contribute to the next section of the dependence graph
			//2) suppose a logical qubit is mapped in mapping 2 but not in mapping 1
				//then that logical qubit must be involved in an adjacent cnot gate
				//therefore, we try mapping common qubits first, and as we consider swaps with unmapped neighbors, we assign qubits
			//3) if a logical qubit is mapped in 1 and 2
				//then we just consider possible swaps of distance 1 away that are close
		 for (auto j = transformedMapping.begin();j != transformedMapping.end();j++) {
			 if (secondMap.find(j->first) == secondMap.end()) {
				 //then this qubit does not contribute to the next maximal mapping 
				 //printf("LOGICAL QUBIT %d DOES NOT CONTRIBUTE TO NEXT MAPPING!\n", j->first);

				 continue;
			 }
			 if (transformedMapping[j->first] != secondMap[j->first]) {
				 pair<int, int> bestSwap = make_pair(-1, -1);
				 int distance = (distances.find(make_pair(j->second, secondMap[j->first])))->second;
				 for (auto i = transformedMapping.begin();i != transformedMapping.end();i++) {
					 if (architectureEdges[j->second].find(i->second) != architectureEdges[j->second].end()) {
						 //then this is a possible swap
						 int firstDistance = (distances.find(make_pair(i->second, secondMap[j->first])))->second;
						 int secondDistance = (distances.find(make_pair(j->second, secondMap[j->first])))->second;

						 if (firstDistance < distance) {
							 //updating swap and distance
							 bestSwap = make_pair(i->first, j->first);
							 distance = firstDistance;
						 }
					 }
				 }

				 //doing the best Swap
				 vector<pair<int, int>> newSwaps = swapsTaken;
				 newSwaps.push_back(bestSwap);
				 //printf("CONSIDERING SWAPPING %d and %d\n", i->first, j->first);

				 //updating transformed mapping
				 map<int, int> updatedMapping = transformedMapping;
				 int tmp = updatedMapping[bestSwap.first];
				 updatedMapping[bestSwap.first] = updatedMapping[bestSwap.second];
				 updatedMapping[bestSwap.second] = tmp;

				 //updating cost
				 int cost = newSwaps.size() + partialMappingStitchingCost(updatedMapping, secondMap, distances, architectureEdges);

				 //enqueueing
				 swapQueue.push(PriorityData(cost, newSwaps, updatedMapping));
			 }
		 }
	 }
	 return results;
 
 }
 


 //this is a swap optimal version of betterLazySwapCircuitBuilder
 void optimalLazySwapCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates) {
	//using priority_queue to push down pathways
	 map<pair<int, int>, int> distances = getDistances(architectureEdges);
	 priority_queue<SecondPriorityData> optimalQueue;

	 //firstly enqueueing the basic arguments-> just the start set
	 optimalQueue.push(SecondPriorityData(0, vector<vector<pair<int, int>>>(), map<int, int>(), map<int, int>(), map<int, int>(), startGates, vector<queue<GateNode*>>()));
	 //below things are for when we break out of the search
	 vector<vector<pair<int, int>>> optimalSwaps;
	 map<int, int> initialMapping;
	 vector<queue<GateNode*>> printQueue;
	 while (optimalQueue.size() > 0) {
		 SecondPriorityData pathway = optimalQueue.top();
		 optimalQueue.pop();
		 

		 //we have 2 special cases -> first map is not set, second map is not set
		 if (pathway.initialMapping.size() == 0) {
			//need to find first mappings, update print queue
			 vector<tuple<map<int, int>, set<int>, set<GateNode*>>> mappings = maximalMapper(architectureEdges, pathway.firstGates);
			 
			 //considering every maximal mapping we get
				//should check if we got a maximal mapping that is perfect or not

		     bool perfectMappingBreak = false;
			 for (auto j = mappings.begin();j != mappings.end();j++) {
				//building the first print queue
				 queue<GateNode*> satisfied= getSatisfiedGates(pathway.firstGates, get<2>(*j));
				 
				 //we have an initial mapping and transformed mapping now
				 //we generate all filled mappings and remove repeats
			 	 vector<queue<GateNode*>> toAdjust = pathway.endingSections;
				 toAdjust.push_back(satisfied);
				 vector<map<int, int>> fullMappings= fillPartialMapping(get<0>(*j), architectureEdges);
				 for (auto z = fullMappings.begin();z != fullMappings.end();z++) {
					 map<int, int> newInitialMapping = *z;	
					 if (get<2>(*j).size() == 0) {
						//we found a perfect mapping, lets just set what we need to and break
						 perfectMappingBreak = true;
						 initialMapping = newInitialMapping;
						 printQueue = toAdjust;
						 break;
					 }
					 optimalQueue.push(SecondPriorityData(0, pathway.swapsTaken, newInitialMapping, pathway.targetMapping, newInitialMapping, get<2>(*j), toAdjust));
				 }
				 if (perfectMappingBreak) {
					 break;
				 }
				// map<int, int> newInitialMapping = get<0>(*j);
				 //we enqueue all these initial pathways with cost 0, because we do not have a swap cost from a single mapping
				 
				 //enqueueing
				 //optimalQueue.push(SecondPriorityData(0, pathway.swapsTaken, newInitialMapping, pathway.targetMapping, newInitialMapping, get<2>(*j), toAdjust));
			 }
			 if (perfectMappingBreak) {
				 break;
			 }
			 continue;
			 
		 } 

		 if (pathway.targetMapping.size() == 0) {
			 //need to find deeper level of maximal mappings
			 //getting cropped architecture
			 map<int, set<int>> croppedArchitecture = getCroppedArchitecture(pathway.initialMapping, architectureEdges);
			 map<pair<int, int>, int> modifiedDistances = getDistances(croppedArchitecture);
			 //now we can start updating the swap cost
			 //printf("ABOUT TO RUN MAXIMAL MAPPER!\n");
			 vector<tuple<map<int, int>, set<int>, set<GateNode*>>> mappings = maximalMapper(croppedArchitecture, pathway.firstGates);
			 for (auto j = mappings.begin();j != mappings.end();j++) {
				 //building the print queue for this step
				 queue<GateNode*> satisfied = getSatisfiedGates(pathway.firstGates, get<2>(*j));
				 //updating the target mapping now
				 map<int, int> newTargetMapping = get<0>(*j);
				 vector<queue<GateNode*>> toAdjust = pathway.endingSections;
				 toAdjust.push_back(satisfied);
				 //getting cost
					//cost only involves heuristic right now since we didnt do any swaps yet
				 int cost = partialMappingStitchingCost(pathway.transformedMapping, newTargetMapping, modifiedDistances, architectureEdges);

				 //enqueuing
				 optimalQueue.push(SecondPriorityData(cost, pathway.swapsTaken, pathway.initialMapping, newTargetMapping, pathway.transformedMapping, get<2>(*j), toAdjust));
			 }
			 continue;
		 }
		 //printf("GOT HERE!\n");
		 //then we stitch the transformed mapping and target mapping together
		 map<int, set<int>> croppedArchitecture = getCroppedArchitecture(pathway.initialMapping, architectureEdges);
		 map<pair<int, int>, int> newDistances = getDistances(croppedArchitecture);
			//if the next gates to start from is empty, we are done
			//we found the optimal swap solution and we can get to printing
			//NOTE: THERE MAY BE MULTIPLE STITCHINGS BETWEEN EACH MAXIMAL PARTITION THAT LEAD TO DIFFERENT RESULT MAPPINGS, WE HAVE TO FIND THE BEST
		 //pair<map<int,int>,vector<pair<int, int>>> newSwaps = stitchMappings(make_pair(pathway.transformedMapping, set<int>()), make_pair(pathway.targetMapping, set<int>()), architectureEdges,distances);
		 //printf("STARTING SWAPS!\n");
		 vector<vector<pair<int, int>>> newSwapList = getSwapPathways(pathway.transformedMapping, pathway.targetMapping, newDistances, croppedArchitecture);
		 //adding the new swaps to the vector of layer swaps
		 bool totalbreak = false;
		 for (auto k = newSwapList.begin();k != newSwapList.end();k++) {
				 //need to apply each set of swaps to the transformed mapping
				 map<int, int> newTransformedMapping = pathway.transformedMapping;
				 for (auto p = (*k).begin();p != (*k).end();p++) {
					 int tmp = newTransformedMapping[p->first];
					 newTransformedMapping[p->first] = newTransformedMapping[p->second];
					 newTransformedMapping[p->second] = tmp;
				 }

				 //enqueue this to priority queue
				 vector<vector<pair<int, int>>> totalSwaps = pathway.swapsTaken;
				 totalSwaps.push_back(*k);
				 //cost will be zero because we want the next iteration to maximal map
				 if (pathway.firstGates.size() == 0) {
					 //then we are done! we found our swap pathway
					 printQueue = pathway.endingSections;
					 optimalSwaps = totalSwaps;
					 initialMapping = pathway.initialMapping;
					 //breaking from loop for printing
					 totalbreak = true;
					 break;
				 }
				 else {
					 
					 optimalQueue.push(SecondPriorityData(0, totalSwaps, pathway.initialMapping, map<int, int>(), newTransformedMapping, pathway.firstGates,pathway.endingSections));
				 }
			 
			 
			 
		 }
		 if (totalbreak) {
			 break;
		 }
		 
	 }
	 map<int, set<int>> croppedArchitecture = getCroppedArchitecture(initialMapping, architectureEdges);
	 scheduleCircuit(printQueue, optimalSwaps, initialMapping, 1, 1, 1, croppedArchitecture);

 }

 //my faster version of the optimal solver
 void betterLazySwapCircuitBuilder(map<int,set<int>> originalArchitectureEdges, set<GateNode*> startGates) {
 //since our initial mapping varies until everything is mapped, we have to print out gates at the end
 //but the process is this: find the best initial maximal mapping (maximal Mapping cost)
	 //gates satisfied for every step of the maximal mapping
	 //distance matrix to be updated at each step
	 map<pair<int, int>, int> distances = getDistances(originalArchitectureEdges);
	// vector<vector<int>> distanceMatrix = getPhysicalDistancesFloydWarshall(originalArchitectureEdges);
	 vector<queue<GateNode*>> maximalMappingSatisfiedGates;
	 //swaps that need to be taken after every maximal partition, if any
	 vector<vector<pair<int, int>>> levelSwaps;
	 //maintaining our initial mapping (every qubit may not be mapped at the start)
	 map<int, int> initialMapping;
	 //lastMapping is tracked in case we need to perform swaps
	 pair<map<int, int>,set<int>> lastMapping;
	 set<GateNode*> currentGates = startGates;
	 map<int, set<int>> architectureEdges = originalArchitectureEdges;
	 while (currentGates.size() > 0) {
		 //doing maximal mapping first
		 printf("STARTING MAXIMAL MAPPING!\n");
		 vector<tuple<map<int,int>,set<int>,set<GateNode*>>> partialMapping=maximalMapper(architectureEdges, currentGates);

		

		 
		 //if this is first maximal mapping we will need to find the best one for remaining cost and set that as our initial mapping
		 if (lastMapping.first.size() == 0) {
			 //finding the best mapping out of these
			 set<GateNode*> endGates;
			 int lowestCost = -1;
			 for (auto z = partialMapping.begin();z != partialMapping.end();z++) {
				 int cost = partialMappingCost(get<0>(*z), get<2>(*z), distances, architectureEdges);
				if (lowestCost == -1 || cost < lowestCost) {
					lastMapping = make_pair(get<0>(*z), get<1>(*z));
					lowestCost = cost;
					endGates = get<2>(*z);
				} 
			 }
				
			 //getting the best filled mapping
			 map<int, int> filledMapping = pickBestFilledMapping(lastMapping.first, architectureEdges, endGates,distances);
			 //printf("GOT OUR BEST FILLED MAPPING!\n");

			 //updating the last mapping
			 set<int> toContinue = lastMapping.second;
			 lastMapping = make_pair(filledMapping, toContinue);

			 //setting our initial mapping to last mapping
			 initialMapping = lastMapping.first;

			 //updating architecture edges
			 architectureEdges = getCroppedArchitecture(filledMapping, architectureEdges);
			 //printf("GOT OUR CROPPED ARCHITECTURE!\n");

			 //updating distance matrix for our new architecture
			 //distanceMatrix = getPhysicalDistancesFloydWarshall(architectureEdges);
			 distances = getDistances(architectureEdges);
			 //printf("UPDATED DISTANCE MATRIX!\n");

			 //adding gates to the print queue
			 queue<GateNode*> toPrint = getSatisfiedGates(currentGates, endGates);

			 //updating currentGates
			 currentGates = endGates;

			//pushing print queue to vector
			maximalMappingSatisfiedGates.push_back(toPrint);
			//printf("FINISHED FIRST ITERATION!\n");
			continue;
		 }
		 else {
			 //printf("STARTED OTHER ITERATION!\n");
			//need to find the closest mapping to the previous mapping, and then consider swaps
			 pair<map<int, int>, set<int>> bestMapping;
			 int bestMappingCost = -1;
			 set<GateNode*> endGates;
			 for (auto z = partialMapping.begin();z != partialMapping.end();z++) {
				 int cost = partialMappingStitchingCost(lastMapping.first, get<0>(*z), distances, architectureEdges);
				 if (bestMappingCost==-1 || cost<bestMappingCost) {
					 //then we can just set this as the best
					 bestMapping = make_pair(get<0>(*z), get<1>(*z));
					 bestMappingCost = cost;
					 endGates = get<2>(*z);
				}
			 }
			 //printf("STARTED SWAPS!\n");
			 //now we find the list of swaps to take
			 pair<map<int, int>, vector<pair<int, int>>> swaps = stitchMappings(lastMapping, bestMapping, architectureEdges,distances);
			 //printf("FINISHED SWAPS\n");

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
			//printf("FINISHED NEXT ITERATION\n");
		 }



	 }
	 
	 //we call schedule here
	
	 scheduleCircuit(maximalMappingSatisfiedGates,levelSwaps,initialMapping,1,1,1,architectureEdges);
	 //printf("REACHED HERE!\n");

	 //if there are unmapped qubits in the initial mapping, then they do not matter and we can assign them whatever we want here
		//specifically, we can assign them any unmapped physical qubit
	 /*
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
	 */
	 /*
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
	 }*/

	 //printing the initial mapping
	 cout << "INITIAL MAPPING!: \n";
	 for (auto z = initialMapping.begin();z != initialMapping.end();z++) {
		 cout << "Logical qubit q" << z->first << " mapped to physical qubit Q" << z->second <<"\n";
	 }
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

 
