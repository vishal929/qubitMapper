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




//this is a big part of my program. it will get partial mappings and the maximal circuit partitions they correspond to
vector <tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>> > maximalMapper( map<int, set<int>> architectureEdges, set<GateNode*> startSet,bool fast,set<GateNode*> totalCompleted);

//this is a helper function to tell us we cropped too many nodes in the maximalMapper, which means at the end of the level, we will get a maximal partition
bool wentTooFar(map<int,set<int>> circuitEdges, map<int,set<int>>architectureEdges,int circuitMaxDegree, int architectureMaxDegree);

//rolls back any changes we made while pruning the architecture in the maximal mapper
map<int,set<int>> rollBack(queue<pair<int,int>>changes, map<int,set<int>>architectureEdges);

//adds a gate edge to the counterpart circuit graph, which will help us prune the architecture graph
map<int,set<int>> addGateEdge(GateNode* chosen, map<int,set<int>>circuitEdges);

//does the actual pruning for finding a mapping for the maximalMapper
 pair<queue<pair<int,int>>,map<int,set<int>>> pruneArchitectureEdges(map<int, set<int>>architectureEdges, unsigned int circuitMinDegree);

 //implementation of a perfect mapper using dfs
 map<int, int> perfectMapper(map<int, set<int>> architectureEdges, set<GateNode*> startSet);

 //gets all maximalMappings of a circuit
 queue<map<int, int>> getMaximalMappings(set<GateNode*> startGates, map<int, set<int>> architectureEdges);

//calculates the maximum distance between any 2 qubits that are mapped in both mappings
 int partialMappingStitchingCost(map<int, int> firstMap, map<int, int> secondMap, map<pair<int,int>,int> distances, map<int, set<int>> architectureEdges);

 //helper for deciding on a maximal partition to fill and pursue initially by deciding an estimated number of swaps based on the gates remaining
 int partialMappingCost(map<int, int> mapping, set<GateNode*> remainingGates, map<pair<int,int>,int> distances, map<int, set<int>> architectureEdges);

 //Uses A* with partialMappingStitchingCost to find an optimal set of swaps to turn a full mapping into the second partial mapping
	//this makes some optimizations for large benchmarks, just so they do not bog
 pair<map<int, int>, vector<pair<int, int>>> stitchMappings(pair<map<int, int>, set<int>> mapOne, pair<map<int, int>, set<int>> mapTwo, map<int, set<int>> architectureEdges,map<pair<int,int>,int> distances);

 //this is my fast mapper (not guarunteed to be optimal)
void betterLazySwapCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates);

//this is a failed optimal mapper
void optimalLazySwapCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates);

//based on a filled mapping, it crops the architecture to only depend on those filled qubits
	//i.e the new architecture will be an induced subgraph of the old architecture, but only with the mapped physical qubits
map<int, set<int>> getCroppedArchitecture(map<int, int> filledMapping, map<int, set<int>> architectureEdges);

//takes a partial mapping and gets all possible fillings for the unmapped qubits
vector<map<int, int>> fillPartialMapping(map<int, int> mapping, map<int, set<int>> architectureEdges);

//considers all possible fillings, and then uses partialMappingCost to decide the best one
map<int, int> pickBestFilledMapping(map<int, int> mapping, map<int, set<int>> architectureEdges, set<GateNode*> remainingGates, map<pair<int,int>,int> distances);

//gets the distances between every 2 physical qubits in the architecture provided as a parameter
map<pair<int, int>, int> getDistances(map<int, set<int>> architectureEdges);
int getNaiveDistance(map<int, set<int>> architectureEdges, int Q1, int Q2);

//this is the pure optimal version of the stitchMappings function
	// this will make no optimizations and considers all swap pathways using A* to stitch our mappings together
	//there might be more than 1 set of swaps that accomplishes this, especially if both mappings are partial, so this is necessary to get an optimal solution
vector<vector<pair<int, int>>> getSwapPathways(map<int, int> one, map<int, int> two, map<pair<int, int>, int> distances, map<int, set<int>> architectureEdges);

//given a level in the given dependence structure, it gets the next level (the result gates do not depend on each other)
set<GateNode*> getNextLevel(set<GateNode*> previousLevel);

//fast mapper, but the circuit prints as mapping occurs, which helps saves up on total space compared to betterLazySwap
void parallelSchedulingFastCircuitBuilder(map<int, set<int>> originalArchitectureEdges, set<GateNode*> startGates);

//FIXED IN HERE!
void improvedParallelSchedulingFastCircuitBuilder(map<int, set<int>> originalArchitectureEdges, set<GateNode*> startGates);


//this is my optimal mapper. every maximal partition is pursued, every swap is considered with A* as the driver to get the best pathways, we pick the pathways with the least # of swaps
void optimalCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates);

queue<GateNode*> getSatisfiedGates(set<GateNode*> beginGates, set<GateNode*> endGates);

void improvedOptimalCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates);


 //priority queue data structure for swaps
	//our cost function is numSwapsSoFar + partialMappingStitchingCost(currentMap, targetMap)
 struct PriorityData {
	 int cost;
	 vector<pair<int, int>> swapsTaken;
	 map<int, int> transformedMapping;

	 PriorityData(int cost, vector<pair<int, int>> swapsTaken, map<int, int> transformedMapping) {
		 this->cost = cost;
		 this->swapsTaken = swapsTaken;
		 this->transformedMapping = transformedMapping;
	 }
	
	 //reversing this operator because we want a min queue
	 bool operator<(const struct PriorityData& toCompare) const{
		 return cost > toCompare.cost; 
	 }

 };

 struct BetterPriorityData {
	 int cost;
	 vector<vector<pair<int, int>>> swapsTaken;
	 map<int, int> transformedMapping;
	 map<int, int> targetMapping;
	 set<GateNode*> firstGates;
	 vector<set<GateNode*>> endingSections;

	 BetterPriorityData(int cost, vector<vector<pair<int,int>>> swapsTaken, map<int,int> transformedMapping, map<int,int> targetMapping, set<GateNode*> firstGates, vector<set<GateNode*>> endingSections) {
		 this->cost = cost;
		 this->swapsTaken = swapsTaken;
		 this->transformedMapping = transformedMapping;
		 this->firstGates = firstGates;
		 this->endingSections = endingSections;
		 this->targetMapping = targetMapping;
	 }

	//reversing this operator because we want a min queue
		 bool operator<(const struct BetterPriorityData& toCompare) const{
			 return cost > toCompare.cost; 
		 }
 };

 //priority queue struct format for optimal swap pathway finder
		//this is specifically for pushing down options in my optimal swap function
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
   //global to keep logical qubit number for ease
   int numLogical = 0;

  //we keep the latencies as globals, just for ease
  int latencySingle;
  int latencyDouble;
  int latencySwap;


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
	//setting global latencies
	latencySingle = latency1;
	latencyDouble = latency2;
	latencySwap = latencySwp;
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
	

	//creating a multimap for quick lookup of architecture edges
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
	
	
	/*Test for lazy circuit builder*/
	//betterLazySwapCircuitBuilder(architectureEdges, firstGates);
	/*optimal test*/
	
	cout << "Please enter 0 if you want the fast implementation (non-optimal), and enter 1 if you want the swap optimal approach\n";
	int x;
	cin >> x;
	if (x == 0) {

		improvedParallelSchedulingFastCircuitBuilder(architectureEdges, firstGates);
	}
	else {

		//optimalCircuitBuilder(architectureEdges, firstGates);
		//optimalLazySwapCircuitBuilder(architectureEdges, firstGates);

		improvedOptimalCircuitBuilder(architectureEdges, firstGates);
	}

	//testing maximal mapper
	/*
	set<GateNode*> currentGates = firstGates;
	map<int, set<int>> architectureToUse = architectureEdges;
	set<GateNode*> encounteredGates;
	map<pair<int, int>, int> distances = getDistances(architectureEdges);
	set<GateNode*> totalCompleted;
	int count = 0;
	while (true) {
		if (currentGates.size() == 0) {
			printf("DONE MAPPING!\n");
			break;
		}
		printf("PRINTING CURRENT FIRST GATES!\n");
		for (auto i = currentGates.begin();i != currentGates.end();i++) {
			cout << "GATE NAME " << (*i)->name << " control: " << (*i)->control << " target: " << (*i)->target << "\n";
		}
		vector<tuple<map<int, int>, set<int>, set<GateNode*>, set<GateNode*>>> mappings = maximalMapper(architectureToUse, currentGates, true,totalCompleted);
		if (count == 0) {
			
			map<int, int> chosenMapping = get<0>(mappings[0]);
			map<int, int> filled = pickBestFilledMapping(chosenMapping, architectureEdges, get<2>(mappings[0]), distances);
			architectureToUse = getCroppedArchitecture(filled, architectureToUse);
		}
		count++;
		printf("PRINTING SATISFIED GATES!\n");
		queue<GateNode*> satisfied;
		for (auto i = get<3>(mappings[0]).begin();i != get<3>(mappings[0]).end();i++) {
			satisfied.push(*i);
			totalCompleted.insert(*i);
		}
		while (satisfied.size() > 0) {
			GateNode* sat = satisfied.front();
			satisfied.pop();
			cout << "GATE NAME " << (sat)->name << " control: " << (sat)->control << " target: " << (sat)->target << "\n";	
		}



		currentGates = get<2>(mappings[0]);
		
		for (auto z = currentGates.begin();z != currentGates.end();) {
			if (totalCompleted.find(*z) != totalCompleted.end()) {
				z = currentGates.erase(z);
			}
			else {
				z++;
			}
		
	}*/
	
	

	//better fast test
	
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
		for (auto i = start.begin();i!=start.end();) {
			
			if (start.find((*i)->controlParent) == start.end() && start.find((*i)->targetParent) == start.end()) {
				//then we add it into queue to get a level order
				startGates.push((*i));
				i=start.erase(i);
			} else{
				i++;
			}
		}
	}

	//going gate by gate through each level, adding edges and then checking if we went too far
	//big note here: gates in each level do not depend on each other, so we know we can always map level by level
	while (startGates.size() > 0) {
		//adding the current level of the dependence graph to our circuit graph and updating node degrees

		GateNode* chosen = startGates.front();
		startGates.pop();
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
		pair<queue<pair<int,int>>,map<int,set<int>>> pruneRes = pruneArchitectureEdges(copyArchitectureEdges,circuitMinDegree);
		queue<pair<int,int>> changes = pruneRes.first;
		copyArchitectureEdges = pruneRes.second;

		unsigned int archMinDegree=copyArchitectureEdges.size();
		for (auto i =copyArchitectureEdges.begin() ; i!=copyArchitectureEdges.end();i++){
			if ((i->second).size() < archMinDegree){
				archMinDegree = (i->second).size();
			}
		}
		
		
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

//just a checker to see if the parent was satisfied in the last maximal mapping or not
bool parentSatisfied(set<GateNode*> startSet, GateNode* parent, GateNode* current) {
	if (current == parent->controlChild) {
		//then we check the targetChild of the parent
		//if we are able to reach start set from the target child, then return true
	}
	else {
		//then we check the controlChild of the parent
		//if we are able to reach a gate in the start set from the control child, then return true
	}
}

//need to see if given a gate and a set of gates, whether the set of gates contains the grandparents of the given gate (this helps for maximal mapper)
bool isGrandParent(GateNode* grandChild, set<GateNode*> startGates, set<GateNode*> gatesToCheck) {
	//go up from the grandchild until we hit a gate in gatesToCheck
	//if we reach a gate in start Gates, then we can stop the upward search
	set<GateNode*> goUpwards;
	set<GateNode*> newSet;
	//first inserting parents of grandChild to the set
	if (grandChild->targetParent != NULL) {
		goUpwards.insert(grandChild->targetParent);
	}
	if (grandChild->controlParent != NULL) {
		goUpwards.insert(grandChild->controlParent);
	}
	if (goUpwards.size() == 0) {
		//printf("START GATE!\n");
		return false;
	}
	//goUpwards.insert(grandChild);
	while (true) {
		for (auto i = goUpwards.begin();i != goUpwards.end();i++) {
			if (gatesToCheck.find(*i) != gatesToCheck.end()) {
				//grandparent confirmed
				return true;
			}
			
		}

		for (auto i = goUpwards.begin();i != goUpwards.end();i++) {
			if (startGates.find(*i) != startGates.end()) {
				//then we did not find grandparent in between
				return false;
			}

			//adding parents to newSet
			if ((*i)->controlParent != NULL) {
				newSet.insert((*i)->controlParent);
			}
			
			if ((*i)->targetParent != NULL) {
				newSet.insert((*i)->targetParent);
			}
		}

		goUpwards = newSet;
		newSet = set<GateNode*>();
	}
	return false;
	


}


//maximal mapper that returns matchings for the next maximal portion of the dependence graph
	//returns a pair
		//the first element is a list of all maximal matchings
		//the second element is jsut a set of mapped architecture nodes --> for quick lookup
		//the third element is the queue of gatenodes which could not be matched and should be considered the "start" for the next iteration
		//this is returned so we can easily start the next iteration of mapping
		//fourth element is the associated set of completed gates
	//the last argument is boolean "fast". If fast is enabled, some speed optimizations are made for compatibility with the larger circuits
	//this may result in a maximal mapping that is not actually maximal, but it is faster for larger circuits
	
vector <tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>> > maximalMapper( map<int, set<int>> architectureEdges, set<GateNode*> startSet,bool fast,set<GateNode*> totalCompleted) {
	//we first grab a copy of the coupling graph and associated nodeDegree vector
	
	map<int,set<int>> copyArchitectureEdges =architectureEdges;
	
	map<int,set<int>> circuitEdges;


	//keeping track of considered gates (we do not want to consider gates twice)
	set<GateNode*> closedSet;
	//this is a queue to insert accounted gates into, will come into play when we actually start mapping
	deque<GateNode*> mappingQueue;
	//keeping track of gates whose requirements have been satisfied
	set<GateNode*> mappingSet;
	set<GateNode*> couldntMap;

		
	queue<GateNode*> startGates;
	set<GateNode*> start = startSet;
		for (auto i = start.begin();i!=start.end();i++) {
			if (start.find((*i)->controlParent) == start.end() && start.find((*i)->targetParent) == start.end() && totalCompleted.find(*i)==totalCompleted.end()) {
				//we only keep the parent
				startGates.push((*i));
			}
			
			
			
		}
	

	//going gate by gate through each level, adding edges and then checking if we went too far
	//big note here: gates in each level do not depend on each other, so we know we can always map level by level
	while (startGates.size() > 0) {
		//adding the current level of the dependence graph to our circuit graph and updating node degrees

		GateNode* chosen = startGates.front();
		startGates.pop();

		if (closedSet.find(chosen) != closedSet.end() || couldntMap.find(chosen) != couldntMap.end()) {
			//no point continuing with this, as we have tried it before
			continue;
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
			}
			if (chosen->targetChild!=NULL){
				GateNode* child = chosen->targetChild;
				if (child->control == chosen->target){
					if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end() || totalCompleted.find(child->targetParent)!=totalCompleted.end()){
						startGates.push(child);
					}
				} else{
					if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end() || totalCompleted.find(child->controlParent)!=totalCompleted.end()){
						startGates.push(child);
					}
				}
				
			}
			
			
			continue;
		}

		//going to closed set and seeing if there are any target, control repeats
			//if so, they are automatically satisfied
		/*
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

				if (chosen->targetChild!=NULL && chosen->targetChild != chosen->controlChild){
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
				repeat = true;
				break;
			}
			
		}
		if (repeat) {
			continue;
		}*/
		
		circuitEdges = addGateEdge(chosen, circuitEdges);
		//cout << "ADDED gate with control " << chosen->control << "target " << chosen->target << "\n";
		//cout << "PRINTING CURRENT CIRCUIT EDGES!\n";
		/*
		for (auto i = circuitEdges.begin();i != circuitEdges.end();i++) {
			cout << i->first << " edges !\n";
			for (auto j = (i->second).begin();j != (i->second).end();j++) {
				cout << *j << ", ";
			}
			cout << "\n";
		}*/
		unsigned int circuitMinDegree = circuitEdges.size();
		
		//unsigned int circuitMinDegree=circuitEdges.size();
		unsigned int circuitMaxDegree=0;
		for (auto i=circuitEdges.begin();i!=circuitEdges.end();i++){
			if ((i->second).size()>circuitMaxDegree){
				circuitMaxDegree = (i->second).size();
			}
			if ((i->second).size()<circuitMinDegree){
				circuitMinDegree = (i->second).size();
			}
		}
		//printf(" MIN DEGREE %d and max degree %d \n", circuitMinDegree, circuitMaxDegree);
		pair<queue<pair<int,int>>,map<int,set<int>>> pruneRes = pruneArchitectureEdges(copyArchitectureEdges,circuitMinDegree);
		queue<pair<int,int>> changes = pruneRes.first;
		copyArchitectureEdges = pruneRes.second;

		unsigned int archMinDegree=copyArchitectureEdges.size();
		for (auto i =copyArchitectureEdges.begin() ; i!=copyArchitectureEdges.end();i++){
			if ((i->second).size() < archMinDegree){
				archMinDegree = (i->second).size();
			}
		}
		
		
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
		//printf("ARCHITECTURE MAX DEGREE: %d\n", architectureMaxDegree);
		if (wentTooFar(circuitEdges,copyArchitectureEdges,circuitMaxDegree,architectureMaxDegree)) {
			//printf("WE WENT TOO FAR!\n");
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
		}
		else {
			//then we can update (add children of the block to the queue
			//adding the chosen gate to the mapping queue
			mappingQueue.push_back(chosen);
			//mappingQueue mirrors the closedSet set
			closedSet.insert(chosen);
			if (chosen->controlChild!=NULL){
				GateNode* child = chosen->controlChild;
				if (child->control == chosen->control){
					//then we need to check target parent
					if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end() || totalCompleted.find(child->targetParent)!= totalCompleted.end() ){
						startGates.push(child);
					}
				} else {
					//then we need to check control parent
					if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end() || totalCompleted.find(child->controlParent)!= totalCompleted.end()){
						startGates.push(child);
					}
				}
			}

			if (chosen->targetChild!=NULL && chosen->targetChild != chosen->controlChild){
				GateNode* child = chosen->targetChild;
				
				if (child->control == chosen->target){
					//then we need to check target parent
					if (child->targetParent == NULL || closedSet.find(child->targetParent)!=closedSet.end() || totalCompleted.find(child->targetParent)!=totalCompleted.end()){
						startGates.push(child);
					}
				} else {
					
					//then we need to check control parent
					if (child->controlParent == NULL || closedSet.find(child->controlParent)!=closedSet.end() || totalCompleted.find(child->controlParent)!=totalCompleted.end()){
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
	//printf("MAPPING QUEUE SIZE: %d\n", max);
	set<GateNode*> debugSet;
while (count != max) {
	GateNode* forMapping = mappingQueue.front();
	mappingQueue.pop_front();
	mappingQueue.push_back(forMapping);
	if ((forMapping->controlParent != NULL && debugSet.find(forMapping->controlParent) == debugSet.end()) || (forMapping->targetParent != NULL && debugSet.find(forMapping->targetParent) == debugSet.end())) {
		if ((forMapping->controlParent != NULL && forMapping->controlParent->control == -1) || (forMapping->targetParent != NULL && forMapping->targetParent->control == -1)) {
		}
	}
	debugSet.insert(forMapping);
	count++;
}


//using a queue to go gate by gate and extend our maximal mappings
queue < tuple<map<int, int>, set<int>, set<GateNode*>>> maximalMappingQueue;
//keeping track of where we are in the mapping

//we do the first iteration
GateNode* firstGate = mappingQueue.front();
mappingQueue.pop_front();
int control = firstGate->control;
int target = firstGate->target;
if (control == -1) {
	//then we will just push an empty mapping onto the queue
	map<int, int> emptyMapping;
	set<int> emptyArch;
	set<GateNode*> initial;
	initial.insert(firstGate);
	maximalMappingQueue.push(make_tuple(emptyMapping, emptyArch, initial));

}
else {
	auto iter = copyArchitectureEdges.begin();
	while (iter != copyArchitectureEdges.end()) {
		if ((iter->second).size() >= circuitEdges[control].size()) {
			//then we found a physical qubit that can accomodate the logical qubit "control" in this maximal mapping
			//we need to check the neighbors of this qubit for possible "target" mappings
			for (auto i = (iter->second).begin();i != (iter->second).end();i++) {
				if (copyArchitectureEdges[*i].size() >= circuitEdges[target].size()) {
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
					maximalMappingQueue.push(make_tuple(suggestedInitialMapping, mappedArchitectureQubits, currentlySatisfied));

				}
			}

		}
		iter++;
	}
}
//starting mapping loop
//keeping track of the best mappings at end of every step
unsigned int lastMapped = 1;
while (mappingQueue.size() > 0) {
	GateNode* toSatisfy = mappingQueue.front();
	mappingQueue.pop_front();


	//we go through all mappings in the queue and see if we can propogate them to the next gate
		//we can easily check if parent gates are satisfied in the current mapping or not with the extra set given
	//queue < tuple<map<int, int>, set<int>,set<GateNode*>>> resultMaximalMappingQueue;

	//maybe do fast optimization here
	//we should cut down on mappings to propogate in order to save on execution time
	//this below step is just for the large benchmarks, without this, they balloon and use too much memory 
		//so with the fast flag enabled, this optimization is made
	while (maximalMappingQueue.size() > 60000 && fast) {
		//throwing out any element
		maximalMappingQueue.pop();
	}

	unsigned int internalStep = 0;
	unsigned int size = maximalMappingQueue.size();
	while (size != 0) {
		tuple<map<int, int>, set<int>, set<GateNode*>> queueEntry = maximalMappingQueue.front();
		maximalMappingQueue.pop();
		size--;

		//seeing if we can propagate the mapping to the next gate, and inserting it in the queue
		map<int, int> currMapping = get<0>(queueEntry);
		set<int> mappedArchQubits = get<1>(queueEntry);
		set<GateNode*> satisfiedGates = get<2>(queueEntry);
		//need to throw out mapping if some conditions are not satisfied
		//this optimization is only made for the large circuits
		if (satisfiedGates.size() < lastMapped && fast) {
			//then we can safely throw this out for now because there are better mappings

			continue;
		}


		/*
		if ((totalCompleted.find(toSatisfy->controlParent) == totalCompleted.end() && satisfiedGates.find(toSatisfy->controlParent) == totalCompleted.end()) ||
			(totalCompleted.find(toSatisfy->targetParent) == totalCompleted.end() && satisfiedGates.find(toSatisfy->targetParent) == totalCompleted.end())){

			//parents not satisfied
			maximalMappingQueue.push(queueEntry);
			continue;
		}*/
			
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
		//updating counter of number of mappings in best mapping so far
		lastMapped = internalStep;
	}
	vector <tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>> > results;
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
	unsigned int numGatesMapped = 0;
	int debugCount = 0;
	while (resultQueue.size() > 0) {
		tuple<map<int, int>, set<int>, set<GateNode*>> queueEntry = resultQueue.front();
		resultQueue.pop();
		if (get<2>(queueEntry).size() < mostMatched) {
			continue;
		}
		

		set<GateNode*> completedGates = get<2>(queueEntry);
		if (debugCount == 0) {
			//printf("MAPPED GATES NUMBER! :%d\n", completedGates.size());
		}
		if (completedGates.size() > numGatesMapped) {
			numGatesMapped = completedGates.size();
		}
		set<GateNode*> finalSet;
		
		//going level by level and seeing what is left to map
		/*
		set<GateNode*> toSkim = startSet;
		while (toSkim.size() > 0) {
			for (auto i = toSkim.begin();i != toSkim.end();) {
				if (completedGates.find(*i) == completedGates.end()) {
					finalSet.insert(*i);
					i = toSkim.erase(i);
				}
				else {
					i++;
				}
			}
			toSkim = getNextLevel(toSkim);
		}*/
		//adding children of completedGates to the end set
		for (auto i = completedGates.begin();i != completedGates.end();i++) {
			if ((*i)->controlChild!=NULL && completedGates.find((*i)->controlChild) == completedGates.end()) {
				//then we add it
				finalSet.insert((*i)->controlChild);
			}

			if ((*i)->targetChild != NULL && completedGates.find((*i)->targetChild) == completedGates.end()) {
				finalSet.insert((*i)->targetChild);
			}
		}

		for (auto i = couldntMap.begin();i != couldntMap.end();i++) {
			finalSet.insert(*i);
		}
		set<GateNode*> finalFinalSet;
		//removing grandparents from the set
		/**/
		for (auto i = finalSet.begin();i != finalSet.end();i++) {
			if (!isGrandParent((*i), startSet, finalSet)) {
				//then the set already contains some ancestor of this gate
				finalFinalSet.insert((*i));
			}
			
		}
		
		results.push_back(make_tuple(get<0>(queueEntry), get<1>(queueEntry), finalFinalSet,completedGates));

		debugCount++;
	}
		
	//now our maximalMappingQueue holds all the generated partial mappings, we have to do some cleanup to get the best ones
		//in addition, we need to do some setup for the future initial gates to consider (because we will run this method multiple times)
	//returning the best partial mappings, and the gates for the next iteration are already modified in firstGates
	//printf("RETURNING RESULTS!\n");
	
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
		//cout << "NO MAPPING POSSIBLE! SPOT 1\n";
		return true;
	}

	if (circuitMaxDegree > architectureMaxDegree) {
		//then we cut out too much from the architecture graph
		//cout << "NO MAPPING POSSIBLE ! SPOT 2\n";
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



 //naive distance implementation for now for distance between 2 physical qubits
 int getNaiveDistance(map<int, set<int>> architectureEdges, int Q1, int Q2) {
		//we basically just use bfs and get shortest distance
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



 //function to get a list of swaps needed to transform one mapping to another
	//returns as the first element of the pair the modified initial mapping (additional logical qubits mapped to mapOne, if necessary)
	//second element of the returned pair is the list of swaps to take
	//if fast is true, then optimization on cutting down priority_queue size is taken
 pair<map<int, int>, vector<pair<int, int>>> stitchMappings(pair<map<int, int>, set<int>> mapOne, pair<map<int, int>, set<int>> mapTwo, map<int, set<int>> architectureEdges, map<pair<int, int>, int> distances) {
	 //initializing a priority queue



	  //priority queue of the next swap to take
	  //the queue is the queue (order) of swaps to take
		 //the map is the modified mapping with the swap
	   
	 priority_queue<PriorityData> swapQueue;

	 map<int, int> firstMap = mapOne.first;
	 map<int, int> transformedMap = firstMap;
	 set<int> firstArchMap = mapOne.second;
	 map<int, int> secondMap = mapTwo.first;
	 set<int> secondArchMap = mapTwo.second;


	 //enqueueing the base state
	 swapQueue.push(PriorityData(partialMappingStitchingCost(firstMap, secondMap, distances, architectureEdges), vector<pair<int, int>>(), transformedMap));
	 //int iterCount = 0;
	 while (swapQueue.size() > 0) {

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
		 //getting the lowest cost path from the queue
		 PriorityData result = swapQueue.top();
		 swapQueue.pop();
		 
		 vector<pair<int, int>> swapsTaken = result.swapsTaken;
		 map<int, int> transformedMapping = result.transformedMapping;
		 
		 //int currCost = result.cost;
		 
		 

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
			 return make_pair(firstMap, swapsTaken);
		 }


		//choosing the closer swaps 
		 for (auto j = transformedMapping.begin();j != transformedMapping.end();j++) {
			 if (secondMap.find(j->first) == secondMap.end()) {
				 //then this qubit does not contribute to the next maximal mapping 

				 continue;
			 }
			 if (transformedMapping[j->first] != secondMap[j->first]) {
				 pair<int, int> bestSwap = make_pair(-1, -1);
				 int distance = (distances.find(make_pair(j->second, secondMap[j->first])))->second;
				 for (auto i = transformedMapping.begin();i != transformedMapping.end();i++) {
					 if (architectureEdges[j->second].find(i->second) != architectureEdges[j->second].end()) {
						 //then this is a possible swap
						 //checking if it is redundant
						 if (swapsTaken.size() > 0) {
							 pair<int, int> lastSwap = swapsTaken[swapsTaken.size() - 1];
							 if ((lastSwap.first == i->first && lastSwap.second == j->first) || (lastSwap.second == i->first && lastSwap.first == j->first)) {
								 //swap is redundant
								 continue;
							 }
						 }
						 
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
				 if (bestSwap.first == -1) {
					 //then no best swap found
					 continue;
				 }
				 vector<pair<int, int>> newSwaps = swapsTaken;
				 newSwaps.push_back(bestSwap);

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
					

 //goes through a pathway of swaps and returns a list of sets, where each set contains all the swaps that can occur at once
	//this will help us generate circuits with lower hits to latency due to swaps
 vector<set<pair<int,int>>> parallelizeSwapPathway(vector<pair<int,int>> swapPathway) {
	//idea: for each swap, go through every other swa	
	 //idea: use a set to determine conflicts (since these swaps are in order)
			//basic loop idea: go through the entire pathway, insert any swaps that do not conflict into a set
	 //basic idea: start at a possible swap, any conflicting swap in the pathway must be in its own set, keep repeating 
	 //printf("GOT HERE!\n");
	 vector<set<pair<int, int>>> parallelizedSwaps;
	 //printf("ALLOCATED VECTOR OF SETS!\n");
	 set<pair<int, int>> addingSet;
	 //printf("ERAI!\n");
	 set<int> currentSet;
	 int count = 0;
	 while (true) {
		 //printf("ITERATION %d: \n",count);
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
				 currentSet.insert(swap.first);
				 currentSet.insert(swap.second);
				 
			 }
			if (currentSet.size() >= numLogical-1) {
					 //cant possibly map more 
					 break;
			}
		 }
		 count++;
		 parallelizedSwaps.push_back(addingSet);
		 //resetting the set
		 addingSet.clear();
		 currentSet.clear();
		 //addingSet = set<pair<int, int>>();
		 //resetting the set of qubits used
		 //currentSet = set<int>();

		 
	 }

	 //now parallelizedSwaps has a vector where each set is a collection of swaps that can be executed in parallel
	 return parallelizedSwaps;
 }

 queue<GateNode*> getSatisfiedGates(set<GateNode*> beginGates, set<GateNode*> endGates) {
	 //we work our way up from endGates until we reach beginGates and then we enqueue
	 queue<GateNode*> finalResult;
	 if (endGates.size() == 0) {
		 //then everything from beginGates to the end is satisfied
		 set<GateNode*> result = beginGates;
		 while (true) {
			 int start = result.size();
			 set<GateNode*> newGates;
			for (auto i = result.begin();i != result.end();i++) {
				if ((*i)->targetChild != NULL) {
					newGates.insert((*i)->targetChild);
				}
				if ((*i)->controlChild != NULL) {
					newGates.insert((*i)->controlChild);
				}
			}

			for (auto i = newGates.begin();i != newGates.end();i++) {
				result.insert(*i);
			}
			if (result.size() == start) {
				//then we hit the end
				for (auto j = result.begin();j != result.end();j++) {
					finalResult.push(*j);
				}
				return finalResult;
			}
		 }
		 
	 }
	 else {
		 //then we work our way up from endGates to beginGates and everything is satisfied in this pathway
		 set<GateNode*> parents;
		 //first adding the parents of endGates (since endGates is not satisfied)
		 for (auto j = endGates.begin();j != endGates.end();j++) {
			 if ((*j)->controlParent != NULL) {
				 parents.insert((*j)->controlParent);
			 }
			 if ((*j)->targetParent != NULL) {
				 parents.insert((*j)->targetParent);
			 }
		 }

		 set<GateNode*> satisfied;
		 while (parents.size() > 0) {
			set<GateNode*> newParents;
			 for (auto i = parents.begin();i != parents.end();i++) {
				 //all parents are satisfied
				 satisfied.insert((*i));
				 //if we hit a start gate we do not propogate
				 if (beginGates.find(*i) == beginGates.end()) {
					 //then we propogate
					 if ((*i)->controlParent != NULL) {
						 newParents.insert((*i)->controlParent);
					 }
					 if ((*i)->targetParent != NULL) {
						 newParents.insert((*i)->targetParent);
					 }
				 }
			}
			 //swapping the sets now
			 parents = newParents;
		 }
		 //adding all satisfied gates to queue and returning 
		 for (auto j = satisfied.begin();j != satisfied.end();j++) {
			 finalResult.push(*j);
		 }
		 return finalResult;
	 }
 }
/*
 //helper function that , when given a list of starting gates, and a list of ending gates, returns a queue of gates up to the ending gates point
 queue<GateNode*> getSatisfiedGates(set<GateNode*> beginGates, set<GateNode*> endGates) {
	 set<GateNode*> start = beginGates;
	 set<GateNode*> closedSet;
	 queue<GateNode*> orderQueue;
	 queue<GateNode*> finalQueue;

	 while (start.size() > 0) {
			for (auto i = start.begin();i!=start.end();) {
				
				if (start.find((*i)->controlParent) == start.end() && start.find((*i)->targetParent) == start.end()) {
					//then we add it into queue to get a level order
					orderQueue.push((*i));
					i=start.erase(i);
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
 }*/

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
	 
	 //returning this new architecture
	 return newArchitecture;
 }

 //helper function to actually schedule the transformed circuit given maximally mapped gates, swaps, and an initial mapping and print the results
	//also prints the initial mapping,cycle depth, and number of gates as requested
	//this takes latencies into account and schedules the gates according to the latency
	//interweave is true if swaps should try and be interleaved in a maximal partition
	//if interweave is false, swaps are just tacked on to the end of a maximal partition
 void scheduleCircuit(vector<set<GateNode*>>& gateMappings,vector<vector<pair<int,int>>> levelSwaps,map<int,int> initialMapping,map<int,set<int>>& architectureEdges) {
	//idea, we just schedule level by level, and whatever is not finished after every cycle, we push to the next	
		//this is the trivial scheduling for the non-interleaved case;
	 //printf("latency single %d\n",latencySingle);
	 //printf("latency double %d\n", latencyDouble);
	 int totalGates = 0;
	 for (auto i = gateMappings.begin();i != gateMappings.end();i++) {
		 totalGates += (*i).size();
	 }
	 printf("CONFIRMING TOTAL GATES IN CIRCUIT: %d \n",totalGates);
	 int numCycles = 0;
	 int numGates = 0;
	 int numSwaps = 0;
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
		 //queue<GateNode*> satisfied = gateMappings[i];
		 //printf("GOT HERE!\n");
		 set<GateNode*> toSchedule=gateMappings[i];
		 set<GateNode*> inProgress;
		 set<GateNode*> finished;
		 //printf("WOOOO!\n");
		 /*
		 while (satisfied.size() > 0) {
			 GateNode* toAdd = satisfied.front();
			 satisfied.pop();
			 toSchedule.insert(toAdd);
		 }*/

		 vector<set<pair<int, int>>> parallelizedSwaps;
		 if (i != gateMappings.size() - 1) {
			 //printf("levelSwaps overall size:%d\n", levelSwaps.size());
			// printf("iteration i: %d\n", i);
			 //printf("levelSwaps size%d\n", levelSwaps.size());
			//printf("STARTING TO PARALLELIZE SWAPS! with levelSwaps[i] size \n",levelSwaps[i].size());
			 parallelizedSwaps = parallelizeSwapPathway(levelSwaps[i]);
			 //printf("SUCCESSFULLY PARALLELIZED SWAPS!\n");
		 }
		 
		 //starting the scheduling
		 int numTotalGates = toSchedule.size();
		 int numSwapsAdded = 0;
		 while (finished.size() < numTotalGates || parallelizedSwaps.size()>0) {
			 //printf("finished size %d\n", finished.size());
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
							 //printf("LATENCY SINGLE! %d\n",latencySingle);
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
							 (toSchedule.find(toConsider->controlParent) == toSchedule.end() && inProgress.find(toConsider->controlParent) == inProgress.end()))) {
						 //then both parents are satisfied, lets try to map these qubits
						 if (logicalSchedule[toConsider->target].second == -1 && logicalSchedule[toConsider->control].second == -1) {
							 //then we can schedule this gate for both qubits
							 //printf("LATENCY DOUBLE: %d\n", latencyDouble);
							 logicalSchedule[toConsider->target] = make_pair(toConsider, latencyDouble);
							 logicalSchedule[toConsider->control] = make_pair(toConsider, latencyDouble);
							//removing this gate from schedule list
							 inProgress.insert(toConsider);
							 //removing gate from toSchedule
							 j = toSchedule.erase(j);
							 //incrementing gate number
							 numGates++;

							 //debug--> making sure that this is valid
							 int firstPhys = toTransform[toConsider->control];
							 int secondPhys = toTransform[toConsider->target];
							 if (architectureEdges[firstPhys].find(secondPhys) == architectureEdges[firstPhys].end()) {
								
								 cout << " SOMETHING WENT TERRIBLY WRONG IN SCHEDULING A QUBIT CNOT!\n";
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
			 if (i != gateMappings.size() - 1) {
				 //printf("STARTING SWAP LOGIC!\n");
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
								 GateNode* newSwapGate = new(GateNode);
								 newSwapGate->name = "swp";
								 newSwapGate->target = y->second;
								 newSwapGate->control = y->first;
								 newSwapGate->targetParent = NULL;
								 newSwapGate->controlParent = NULL;

								 //adding this gate to the schedule and removing it from this set
								 logicalSchedule[y->first] = make_pair(newSwapGate, latencySwap);
								 logicalSchedule[y->second] = make_pair(newSwapGate, latencySwap);
								 y = (*z).erase(y);

								 //incrementing number of gates
								 numGates++;
								 //incrementing number of swaps added
								 numSwapsAdded++;

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
						 z = parallelizedSwaps.erase(z);
					 }
				 }
			 }


			 //we have scheduled all that we can right now, let us move onto the decrement step for cycle logic
				//we decrement all in progress gates by 1
			 for (int z = 0; z < numLogical;z++) {
				 //printf("CYCLE ITERATION %d\n", numCycles);
				if (logicalSchedule[z].second > 0) {
					//decrementing
					GateNode* work = logicalSchedule[z].first;
					int cyclesLeft = logicalSchedule[z].second;
					cyclesLeft--;
					 if (cyclesLeft == 0) {
						 //then this gate is done executing, we should remove it from the inProgress set and from logicalSchedule
						 inProgress.erase(work);
						 finished.insert(work);
						 //if this gate is a swap gate, then we transform the mapping
						 //need to make both logical and control qubits the dummy
						 logicalSchedule[work->target] = make_pair(dummy, -1);
						 if (work->control != -1) {
							 logicalSchedule[work->control] = make_pair(dummy, -1);
						 }
						 if (work->name == "swp") {
							 numSwaps++;
							 //we do the swap change
							 //debug --> checking if swap is valid or not
							 int physOne = toTransform[work->control];
							 int physTwo = toTransform[work->target];
							 
							 if (architectureEdges[physOne].find(physTwo) == architectureEdges[physOne].end()) {
								cout << "SOMETHING WENT TERRIBLY WRONG WITH MAPPING FOR SWAP!\n";
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
					 else {
						 logicalSchedule[z] = make_pair(work, cyclesLeft);
					 }
					
				}
			 }
			 //printf("DONE WITH CYCLE!\n");
			 
			 //cycle is over 
			 numCycles++;
		 }
		 
	 }
	 cout << "NUMBER OF GATES: " << numGates << "\n";
	 cout << "NUMBER OF SWAPS: " << numSwaps << "\n";
	 cout << "NUMBER OF CYCLES: " << numCycles<<"\n";

	 cout << "INITIAL MAPPING: \n";
	 for (auto i = initialMapping.begin();i != initialMapping.end();i++) {
		 cout << "Logical q[" << i->first << "] mapped to physical Q[" << i->second << "]\n";
	 }

	 	
	 
 }

 vector<vector<pair<int, int>>> getSwapPathways(map<int,int> one, map<int,int> two, map<pair<int,int>,int> distances,map<int,set<int>> architectureEdges) {
	 //printf("STARTED SWAP PATHWAYS!\n");
	priority_queue<PriorityData> swapQueue;
	vector <vector<pair<int, int>>> results;
	 map<int, int> firstMap = one;
	 map<int, int> transformedMap = one;
	 map<int, int> secondMap = two;
		
	 /*
	 printf("GET SWAP PATHWAYS!\n");
	 printf("FIRST MAP!\n");
	 for (auto i = one.begin();i != one.end();i++) {
		 cout << i->first << ", " << i->second << "\n";
	 }

	 printf("SECOND MAP!\n");
	 for (auto i = two.begin();i != two.end();i++) {
		 cout << i->first << ", " << i->second << "\n";
	 }*/


	 //enqueueing the base state
	 swapQueue.push(PriorityData(partialMappingStitchingCost(firstMap, secondMap, distances, architectureEdges), vector<pair<int, int>>(), transformedMap));
	 //int iterCount = 0;
	 while (swapQueue.size() > 0) {
		 //optimization from stitchMapping commented out, in case it should be needed here
		 /*
		 while (swapQueue.size() > 60000 ) {
			 //just an optimization to cut down on mappings "we think" are bad, especially useful for large benchmarks
			 //we will only keep the top 3 mappings from swapQueue
			 priority_queue<PriorityData> tmp;
			 for (int i = 0;i < 3;i++) {
				 tmp.push(swapQueue.top());
				 swapQueue.pop();
			 }
			 swapQueue = tmp;
		 }*/
		 //getting the lowest cost path from the queue
		 PriorityData result = swapQueue.top();
		 swapQueue.pop();
		 
		 vector<pair<int, int>> swapsTaken = result.swapsTaken;
		 map<int, int> transformedMapping = result.transformedMapping;
		 
		 //int currCost = result.cost;
		 

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
			 //printf("HIT SWAP EQUAL!\n");
			 //then we hit the end
			 //returning the initial mapping, and list of swaps to take
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


		//considering swaps that take the unagreeing qubits between both maps closer	
		 for (auto j = transformedMapping.begin();j != transformedMapping.end();j++) {
			 if (secondMap.find(j->first) == secondMap.end()) {
				 //then this qubit does not contribute to the next maximal mapping 

				 continue;
			 }
			 if (transformedMapping[j->first] != secondMap[j->first]) {
				// pair<int, int> bestSwap = make_pair(-1, -1);
				 //int distance = (distances.find(make_pair(j->second, secondMap[j->first])))->second;
				 for (auto i = transformedMapping.begin();i != transformedMapping.end();i++) {
					 if (architectureEdges[j->second].find(i->second) != architectureEdges[j->second].end()) {
						 //then this is a possible swap
						 //checking if this is redundant
						 if (swapsTaken.size() > 0) {
							 pair<int, int> lastSwap = swapsTaken[swapsTaken.size() - 1];
							 if ((lastSwap.first == j->first && lastSwap.second == i->first) || (lastSwap.second == j->first && lastSwap.first == i->first)) {
								 //redundant swap
								 continue;
							}
						 }
						 //we consider every possible 1 swap here which gets us closer to the end state
						 /*
						 vector<pair<int, int>> newSwaps = swapsTaken;
							 newSwaps.push_back(make_pair(i->first, j->first));
							 map<int, int> updatedMapping = transformedMapping;
							 int tmp = updatedMapping[j->first];
							 updatedMapping[j->first] = updatedMapping[i->first];
							 updatedMapping[i->first] = tmp;
							 int cost = newSwaps.size() + partialMappingStitchingCost(updatedMapping, secondMap, distances, architectureEdges);
							 swapQueue.push(PriorityData(cost, newSwaps, updatedMapping));
							*/
						 int firstDistance = (distances.find(make_pair(i->second, secondMap[j->first])))->second;
						 int secondDistance = (distances.find(make_pair(j->second, secondMap[j->first])))->second;
						
						 if (firstDistance < secondDistance) {
							 //doing every swap which brings us closer
							 //updating swap and distance
							// bestSwap = make_pair(i->first, j->first);
							 //distance = firstDistance;
							 vector<pair<int, int>> newSwaps = swapsTaken;
							 newSwaps.push_back(make_pair(i->first, j->first));
							 map<int, int> updatedMapping = transformedMapping;
							 int tmp = updatedMapping[j->first];
							 updatedMapping[j->first] = updatedMapping[i->first];
							 updatedMapping[i->first] = tmp;
							 int cost = newSwaps.size() + partialMappingStitchingCost(updatedMapping, secondMap, distances, architectureEdges);
							 swapQueue.push(PriorityData(cost, newSwaps, updatedMapping));

						 }
						 
					 }
				 }
				/*
				 //doing the best Swap
				 vector<pair<int, int>> newSwaps = swapsTaken;
				 newSwaps.push_back(bestSwap);

				 //updating transformed mapping
				 map<int, int> updatedMapping = transformedMapping;
				 int tmp = updatedMapping[bestSwap.first];
				 updatedMapping[bestSwap.first] = updatedMapping[bestSwap.second];
				 updatedMapping[bestSwap.second] = tmp;

				 //updating cost
				 int cost = newSwaps.size() + partialMappingStitchingCost(updatedMapping, secondMap, distances, architectureEdges);

				 //enqueueing
				 swapQueue.push(PriorityData(cost, newSwaps, updatedMapping));
				 */
			 }
		 }
	 }
	 return results;
 
 }

 void improvedOptimalCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates) {
	//we first get every single initial mapping from each maximal mapping	
	 //for each initial mapping, we initialize a priority queue with transformed mapping and swaps etc
		//so for each maximal mapping and initial mapping, we are able to minimize the swaps
		//therefore, we exhaust all the initial options and then we are done
	 queue < pair<map<int, int>, priority_queue<BetterPriorityData>>> mappingQueue;
	 vector<vector<pair<int, int>>> bestSwaps;
	 map<int, int> bestInitialMapping;
	 vector<set<GateNode*>> bestEndingSections;
	 //printf("INITIAL SIZE: %d\n", bestEndingSections.size());
	 //just an indicator for the first part

	 //doing initial setup
	 vector<tuple<map<int, int>, set<int>, set<GateNode*>, set<GateNode*>>> results = maximalMapper(architectureEdges, startGates, false, set<GateNode*>());
	 for (auto j = results.begin();j != results.end();j++) {
		 vector<map<int, int>> filledMappings = fillPartialMapping(get<0>(*j), architectureEdges);
		 set<GateNode*> ending = get<2>(*j);
		 set<GateNode*> completed = get<3>(*j);
		 for (auto i = ending.begin();i != ending.end();) {
			 if (completed.find(*i) != completed.end()) {
				 i = ending.erase(i);
			 }
			 else {
				 i++;
			 }
		 }
		 //checking for perfect mapping
			
		 if (ending.size() == 0) {
			 printf("PERFECT MAPPING!\n");
			 vector<set<GateNode*>> finished;
			 finished.push_back(completed);
			 scheduleCircuit(finished, vector<vector<pair<int, int>>>(), get<0>(*j), architectureEdges);
			 return;
		 }
		 //adding every mapping and last gates set to the queue
		 for (auto z = filledMappings.begin();z != filledMappings.end();z++) {
			 map<int, set<int>> croppedArch = getCroppedArchitecture(*z, architectureEdges);
			 map<pair<int, int>, int> distances = getDistances(croppedArch);
			 //getting a second mapping and adding to the queue with heuristic cost
			 vector<tuple<map<int, int>, set<int>, set<GateNode*>, set<GateNode*>>> furtherResults = maximalMapper(croppedArch, ending, false, completed);
			 for (auto e = furtherResults.begin();e != furtherResults.end();e++) {
				 //adding the second mappings to the queue with cost
				 /*
				 printf("GOT SECOND MAPPING!\n");
				 for (auto l = get<0>(*e).begin();l != get<0>(*e).end();l++) {
					 cout << l->first << " , " << l->second << "\n";
				 }*/
				 vector<set<GateNode*>> toAdd;
				 toAdd.push_back(completed);
				 toAdd.push_back(get<3>(*e));
				 /*
				 int numGates = 0;
				 for (auto h = toAdd.begin();h != toAdd.end();h++) {
					 numGates += (*h).size();
				 }
				 printf("INITIAL ENDING GATES PUSHED %d\n", numGates);
				 */
				 set<GateNode*> end = get<2>(*e);
				 for (auto h = end.begin();h != end.end();) {
					 if (completed.find(*h) != completed.end() || get<3>(*e).find(*h) != get<3>(*e).end()) {
						 h = end.erase(h);
					 }
					 else {
						 h++;
					 }
				 }
				
					priority_queue<BetterPriorityData> deeperQueue;
					int heuristic = partialMappingStitchingCost(*z, get<0>(*e), distances, croppedArch);
					map<int, int> initialMap = *z;
					map<int, int> targetMap = get<0>(*e);
					//printf("TARGET MAP SIZE: %d\n", targetMap.size());
					deeperQueue.push(BetterPriorityData(heuristic, vector<vector<pair<int, int>>>(), initialMap,targetMap, end, toAdd));
					mappingQueue.push(make_pair(initialMap, deeperQueue));
				 
			 }
			 
		 }
	 }

	 //printf("INITIAL MAPPING QUEUE SIZE! :%d\n",mappingQueue.size());

	 //setup is over, we will now start the double queue stuff
	 while (mappingQueue.size() > 0) {
		 //printf("ITERATION!\n");
		 //completely exhausting all inner mappings here first
		 pair<map<int, int>, priority_queue<BetterPriorityData>> popped = mappingQueue.front();
		 mappingQueue.pop();
		 map<int, set<int>> croppedArch = getCroppedArchitecture(popped.first, architectureEdges);
		 map<pair<int, int>, int> distances = getDistances(croppedArch);

		 priority_queue<BetterPriorityData> deeperQueue = popped.second;
		 
		 while (deeperQueue.size() > 0) {
			 //exhausting all possibilities for this mapping
			 BetterPriorityData innerPopped = deeperQueue.top();
			 deeperQueue.pop();

			 //printf("INNER POPPED TARGET MAP SIZE: %d\n", innerPopped.targetMapping.size());
			 //printf("INNER POPPED TRANSFORMED MAP SIZE: %d\n", innerPopped.transformedMapping.size());

			 

			 if (innerPopped.firstGates.size() == 0 && innerPopped.targetMapping.size()==0) {
				 //then we are done, we found an optimal swap pathway for this given initial mapping
				 int numSwaps = 0;
				 for (auto z = innerPopped.swapsTaken.begin();z != innerPopped.swapsTaken.end();z++) {
					 numSwaps += (*z).size();
				 }
				 int numBestSwaps = 0;
				 for (auto z = bestSwaps.begin();z != bestSwaps.end();z++) {
					 numBestSwaps += (*z).size();
				 }
				 if (bestEndingSections.size() == 0) {
					 //then we just fill
					 bestSwaps = innerPopped.swapsTaken;
				     bestInitialMapping = popped.first;
					 bestEndingSections = innerPopped.endingSections;
				 }
				 else if (numSwaps < numBestSwaps) {
					 //we found better swap pathway
					 bestSwaps = innerPopped.swapsTaken;
				     bestInitialMapping = popped.first;
					 bestEndingSections = innerPopped.endingSections;
				 }
				 int numGates = 0;
				 for (auto z = innerPopped.endingSections.begin();z != innerPopped.endingSections.end();z++) {
					 numGates += (*z).size();
				 }
				 //printf("FOUND GOOD PATH WITH NUM SWAPS %d and num sections %d, and total gates %d\n", numSwaps,innerPopped.endingSections.size(),numGates);
				 break;
			 }


			 //else we call maximal mapper, stitch the swaps and enqueue
			 set<GateNode*> completedSoFar; 
			 for (auto i = innerPopped.endingSections.begin();i != innerPopped.endingSections.end();i++) {
				 for (auto z = (*i).begin();z != (*i).end();z++) {
					 completedSoFar.insert(*z);
				 }
			 }

			 //then we stitch here and then the target mapping is empty
			 /*
			 printf("INNER SECOND MAPPING!\n");
			 for (auto l = innerPopped.targetMapping.begin();l != innerPopped.targetMapping.end();l++) {
				 cout << l->first << ", " << l->second << "\n";
			 }*/
			 vector<vector<pair<int, int>>> swapPathways = getSwapPathways(innerPopped.transformedMapping, innerPopped.targetMapping, distances, croppedArch);
			 //printf("ENTIRE SWAP PATHWAY SIZE: %d\n", swapPathways.size());
			 //performing swaps to transform our mapping
			 for (auto z = swapPathways.begin();z != swapPathways.end();z++) {
				 //printf("SIZE OF SWAP PATHWAY: %d\n", (*z).size());
				 vector<vector<pair<int, int>>> newSwapsTaken = innerPopped.swapsTaken;
				 newSwapsTaken.push_back(*z);
				 map<int, int> toTransform = innerPopped.transformedMapping;
				 for (auto y = (*z).begin();y != (*z).end();y++) {
					// cout << "SWAP " << y->first << y->second << "\n";
					 int tmp = toTransform[y->first];
					 toTransform[y->first] = toTransform[y->second];
					 toTransform[y->second] = tmp;
				 }
				 //enqueuing from this pathway
				 int numSwaps = 0;
				 for (auto i = newSwapsTaken.begin();i != newSwapsTaken.end();i++) {
					 numSwaps += (*i).size();
				 }
				 //getting the next maximal mapping  if lastGates is not empty
				 if (innerPopped.firstGates.size() != 0) {
					vector<tuple<map<int, int>, set<int>, set<GateNode*>, set<GateNode*>>> results = maximalMapper(croppedArch, innerPopped.firstGates, false, completedSoFar);

					//enqueuing each with heuristic cost + numSwaps
					 for (auto h = results.begin();h != results.end();h++) {
						 set<GateNode*> ending = get<2>(*h);
						 set<GateNode*> completed = get<3>(*h);
						 for (auto z = ending.begin();z != ending.end();) {
							 if (completed.find(*z) != completed.end() || completedSoFar.find(*z) != completedSoFar.end()) {
								 z = ending.erase(z);
							 }
							 else {
								 z++;
							 }
						 }

						 vector<set<GateNode*>> endings = innerPopped.endingSections;
						 endings.push_back(completed);
						 int heuristic = partialMappingStitchingCost(toTransform, get<0>(*h), distances, croppedArch);
						 deeperQueue.push(BetterPriorityData(numSwaps + heuristic, newSwapsTaken, toTransform, get<0>(*h), ending, endings));
					 }
				 }
				 else {
					 //we enqeueu with empty mapping and num swaps
					 deeperQueue.push(BetterPriorityData(numSwaps, newSwapsTaken, toTransform, map<int, int>(), innerPopped.firstGates, innerPopped.endingSections));
				 }
				 
				 
			 }
			 

			 
			 

		 }
	 }
	 //printf("FINAL SECTIONS SIZE %d\n", bestEndingSections.size());

	 //now we have the best mapping pathway possible
	 map<int, set<int>> croppedArch = getCroppedArchitecture(bestInitialMapping, architectureEdges);
	 scheduleCircuit(bestEndingSections, bestSwaps, bestInitialMapping, croppedArch);


	 

 }

/*
 void optimalCircuitBuilder(map<int, set<int>> architectureEdges, set<GateNode*> startGates) {
	//using priority_queue to push down pathways
	 map<pair<int, int>, int> distances = getDistances(architectureEdges);
	 priority_queue<SecondPriorityData> optimalQueue;

	 //firstly enqueueing the basic arguments-> just the start set
	 //optimalQueue.push(SecondPriorityData(0, vector<vector<pair<int, int>>>(), map<int, int>(), map<int, int>(), map<int, int>(), startGates, vector<queue<GateNode*>>()));
	 //below things are for when we break out of the search
	 vector<vector<pair<int, int>>> optimalSwaps;
	 map<int, int> initialMapping;
	 vector<queue<GateNode*>> printQueue;

	 //initializing the priorityQueue with the first 2 sets of mappings (target and initial)
	 vector<tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>>> result = maximalMapper(architectureEdges, startGates, false,set<GateNode*>());
	 //checking for perfect initial mapping
	 if (get<2>(result[0]).size() == 0) {
		 printf("PERFECT MAPPING!\n");
		 //perfect mapping
		 map<int, set<int>> croppedArch = getCroppedArchitecture(get<0>(result[0]), architectureEdges);
		 queue<GateNode*> mapped = getSatisfiedGates(startGates, get<2>(result[0]));
		 printQueue.push_back(mapped);
		 scheduleCircuit(printQueue, optimalSwaps, get<0>(result[0]), croppedArch);
		 return;
	 }

	 //else for each partial mapping, we get fillings, and for each filling we get a second mapping and from there we stitch and repeat
	 for (auto i = result.begin();i != result.end();i++) {
		 vector<map<int, int>> filledMappings = fillPartialMapping(get<0>(*i), architectureEdges);
		 queue<GateNode*> satisfied = getSatisfiedGates(startGates, get<2>(*i));
		 for (auto j = filledMappings.begin();j != filledMappings.end();j++) {
			 //for each filled mapping, we crop the architecture, update distance, and do the next set of mappings
			 map<int,set<int>> croppedArch = getCroppedArchitecture(*j, architectureEdges);
			 vector<tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>>> deepResults = maximalMapper(croppedArch, get<2>(*i), false,set<GateNode*>());
			 for (auto z = deepResults.begin();z != deepResults.end();z++) {
				//now we can get the stitching distance and enqueue into priority queue
				 map<pair<int, int>, int> distances = getDistances(croppedArch);
				 queue<GateNode*> deepSatisfied = getSatisfiedGates(get<2>(*i), get<2>(*z));
				 int cost = partialMappingStitchingCost(*j, get<0>(*z), distances, croppedArch);
				 vector<queue<GateNode*>> bigQueue;
				 bigQueue.push_back(satisfied);
				 bigQueue.push_back(deepSatisfied);
				 optimalQueue.push(SecondPriorityData(cost, vector<vector<pair<int, int>>>(), *j, get<0>(*z), *j, get<2>(*z), bigQueue));
			 }
			
		 }
	 }

	 //finished initializing the priority queue, now lets swap!
	 while (optimalQueue.size() > 0) {
		 SecondPriorityData popped = optimalQueue.top();
		 optimalQueue.pop();
		 map<int, set<int>> croppedArch = getCroppedArchitecture(popped.initialMapping,architectureEdges);
		 map<pair<int, int>, int> distances = getDistances(croppedArch);
		 //if the second mapping is empty, then this is optimal pathway
		 if (popped.targetMapping.size() == 0) {
			 //then we  are done
			 scheduleCircuit(popped.endingSections, popped.swapsTaken, popped.initialMapping, croppedArch);
			 return;
		 }

		 //else we stitch the transformed to the target, update the cost, call maximalMapper again, and continue
		 //stitching now
		 map<int, int> transformedMap = popped.transformedMapping;
		 map<int, int> dest = popped.targetMapping;
		 vector<vector<pair<int, int>>> allSwapPaths = getSwapPathways(transformedMap, dest, distances, croppedArch);
		 for (auto i = allSwapPaths.begin();i != allSwapPaths.end();i++) {
			 //going through each pathway, updating the transformedMap
			 map<int, int> newTransformedMap=popped.transformedMapping;
			 vector<vector<pair<int, int>>> swapPathWay = popped.swapsTaken;
			 swapPathWay.push_back(*i);
			 for (auto z = (*i).begin();z != (*i).end();z++) {
				 //transforming the map
				 int tmp = newTransformedMap[z->first];
				 newTransformedMap[z->first] = newTransformedMap[z->second];
				 newTransformedMap[z->second] = tmp;
			 }
			 //now map is transformed, we need to go through all second mappings, if there are any
			 if (popped.firstGates.size() == 0) {
				 //theres nothing left, we can just enqueue this back with an empty target
				 int totalSwaps = 0;
				 for (auto z = swapPathWay.begin();z != swapPathWay.end();z++) {
					 totalSwaps += (*z).size();
				 }
				 //only our total swaps is the cost
				 optimalQueue.push(SecondPriorityData(totalSwaps, swapPathWay, popped.initialMapping, map<int,int>(), newTransformedMap, popped.firstGates,popped.endingSections));
			 }
			 else {
				 //we proceed with maximal mapping
				 vector<tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>>> res = maximalMapper(croppedArch, popped.firstGates, false,set<GateNode*>());
				 //now we have to go through each one, compute the cost and the satisfied gates and enqueue them
				 for (auto z = res.begin();z != res.end();z++) {
					 queue<GateNode*> satisfied = getSatisfiedGates(popped.firstGates, get<2>(*z));
					 vector<queue<GateNode*>> satisfiedSoFar = popped.endingSections;
					 satisfiedSoFar.push_back(satisfied);
					 int numSwaps = 0;
					 for (auto p = swapPathWay.begin();p != swapPathWay.end();p++) {
						 numSwaps += (*p).size();
					 }
					 int cost = numSwaps + partialMappingStitchingCost(newTransformedMap, get<0>(*z), distances, croppedArch);
					 optimalQueue.push(SecondPriorityData(cost, swapPathWay, popped.initialMapping, get<0>(*z), newTransformedMap, get<2>(*z), satisfiedSoFar));
				 }
			 }
		 }
	 }
 }*/

/*
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
		 if (pathway.firstGates.size()==0) {
			 //we are done here
			 optimalSwaps = pathway.swapsTaken;
			 printQueue = pathway.endingSections;
			 initialMapping = pathway.initialMapping;
			 //printf("REACHED!\n");
			 break;
		 }
		 

		 //we have 2 special cases -> first map is not set, second map is not set
		 if (pathway.initialMapping.size() == 0) {
			//need to find first mappings, update print queue
				//fast is set to false for optimality
			 vector<tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>>> mappings = maximalMapper(architectureEdges, pathway.firstGates, false,set<GateNode*>());
			 
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
					 //start with cost 0 so that all these get propogated and account for swaps
					 optimalQueue.push(SecondPriorityData(0, pathway.swapsTaken, newInitialMapping, pathway.targetMapping, newInitialMapping, get<2>(*j), toAdjust));
				 }
				 if (perfectMappingBreak) {
					 break;
				 }
				 
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
			 vector<tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>>> mappings = maximalMapper(croppedArchitecture, pathway.firstGates,false,set<GateNode*>());
			 for (auto j = mappings.begin();j != mappings.end();j++) {
				 //building the print queue for this step
				 queue<GateNode*> satisfied = getSatisfiedGates(pathway.firstGates, get<2>(*j));
				 //updating the target mapping now
				 map<int, int> newTargetMapping = get<0>(*j);
				 vector<queue<GateNode*>> toAdjust = pathway.endingSections;
				 toAdjust.push_back(satisfied);
				 //getting cost
					//cost only involves heuristic right now since we didnt do any swaps yet
					//cost is zero right now because no stitching involved
				 //cost is however many swaps it took to get here + heuristic
				 vector<vector<pair<int,int>>> allSwaps = pathway.swapsTaken;
				 int numSwaps = 0;
				 for (auto i = allSwaps.begin();i != allSwaps.end();i++) {
					 numSwaps += (*i).size();
				 }
				 int cost = numSwaps +partialMappingStitchingCost(pathway.transformedMapping, newTargetMapping, modifiedDistances, croppedArchitecture) ;

				 //enqueuing
				 optimalQueue.push(SecondPriorityData(cost, pathway.swapsTaken, pathway.initialMapping, newTargetMapping, pathway.transformedMapping, get<2>(*j), toAdjust));
			 }
			 continue;
		 }
		 //then we stitch the transformed mapping and target mapping together
		 map<int, set<int>> croppedArchitecture = getCroppedArchitecture(pathway.initialMapping, architectureEdges);
		 map<pair<int, int>, int> newDistances = getDistances(croppedArchitecture);
			//if the next gates to start from is empty, we are done
			//we found the optimal swap solution and we can get to printing
			//NOTE: THERE MAY BE MULTIPLE STITCHINGS BETWEEN EACH MAXIMAL PARTITION THAT LEAD TO DIFFERENT RESULT MAPPINGS, WE HAVE TO FIND THE BEST
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
				 
					 //need to determine number of swaps taken for cost 
					 int numSwaps = 0;
					 for (auto j = totalSwaps.begin();j != totalSwaps.end();j++) {
						 numSwaps += (*j).size();
					 }
					 //printf("NUM SWAPS %d\n", numSwaps);
					 optimalQueue.push(SecondPriorityData(numSwaps, totalSwaps, pathway.initialMapping, map<int, int>(), newTransformedMapping, pathway.firstGates,pathway.endingSections));
				 
			 
			 
			 
		 }
		 if (totalbreak) {
			 break;
		 }
		 
	 }
	 //printf("REACHED HERE!\n");
	 map<int, set<int>> croppedArchitecture = getCroppedArchitecture(initialMapping, architectureEdges);
	 //printf("AFTER CROP!\n");
	 int numSwaps = 0;
	 for (auto i = optimalSwaps.begin();i != optimalSwaps.end();i++) {
		 numSwaps += (*i).size();
	 }
	 //printf(" WE HAVE %d optimal swaps\n", numSwaps);
	 scheduleCircuit(printQueue, optimalSwaps, initialMapping,  croppedArchitecture);

	 //printing the initial mapping
	 cout << "INITIAL MAPPING: \n";
	 for (auto i = initialMapping.begin();i != initialMapping.end();i++) {
		 cout << "LOGICAL q" << i->first << " mapped to physical Q" << i->second<<"\n";
	 }

 }*/

void improvedParallelSchedulingFastCircuitBuilder(map<int, set<int>> originalArchitectureEdges, set<GateNode*> startGates) {
	//will eventually need to print this out
	 int numSwaps = 0;
	 int numGates = 0;
	 int numCycles = 0;
	 map<pair<int, int>, int> distances = getDistances(originalArchitectureEdges);
	 set<GateNode*> totalCompleted;
	 //queue<GateNode*> forSchedule;
	 set<GateNode*> toSchedule;
	 set<GateNode*> inProgress;
	 set<GateNode*> finished;
	 vector<pair<GateNode*, int>> logicalSchedule;
	 GateNode* dummy = new GateNode();
	 for (int i = 0;i < numLogical;i++) {
		 //-1 means qubit is free to participate in some mapping
		 logicalSchedule.push_back( make_pair(dummy, -1));
	 }
	 //maintaining our initial mapping (every qubit may not be mapped at the start)
	 map<int, int> initialMapping;
	 //keeping track of our current state
	 map<int, int> toTransform;
	 //keeping track of the target mapping for swaps
	 map<int, int> targetMapping;
	 set<GateNode*> currentGates = startGates;
	 map<int, set<int>> architectureEdges = originalArchitectureEdges;

	 while (currentGates.size() > 0) {
		 if (totalCompleted.size() > 1000) {
			 set<GateNode*> otherSet;
			 //optimization just to reduce memory usage

			 //if a gate has both its children in this set, then we can remove it
			 //on the other hand, if a node does not have both children in the set, then we can keep it
			 for (auto i = totalCompleted.begin();i != totalCompleted.end();i++) {
				 bool notBoth = false;
				 if (((*i)->controlChild != NULL && totalCompleted.find((*i)->controlChild) == totalCompleted.end()) || ((*i)->targetChild != NULL && totalCompleted.find((*i)->targetChild) == totalCompleted.end())) {
					 otherSet.insert(*i);
				 }
				
			 }

			 //swapping sets
			 totalCompleted = otherSet;
		 }

		 if (initialMapping.size() == 0) {
			 //need to get initial maximal mapping
			 vector<tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>>> partialMappings = maximalMapper(architectureEdges, currentGates, true,totalCompleted);
			 //picking lowest cost from these
			 map<int, int> bestMapping;
			 set<GateNode*> endGates;
			 int cost = -1;
			 for (auto i = partialMappings.begin();i != partialMappings.end();i++) {
				 int mapCost = partialMappingCost(get<0>(*i), get<2>(*i), distances, architectureEdges);
				 if (cost == -1 || mapCost < cost) {
					 bestMapping = get<0>(*i);
					 cost = mapCost;
					 endGates = get<2>(*i);
					 toSchedule = get<3>(*i);
				 }
			 }

			 //filling the initial mapping
			 map<int, int> fullInitialMapping = pickBestFilledMapping(bestMapping, architectureEdges, endGates, distances);
			 initialMapping = fullInitialMapping;
			 
			 toTransform = fullInitialMapping;
			 //forSchedule = getSatisfiedGates(currentGates, endGates);
			 currentGates = endGates;
			 

			 //adding to completed gates
			 for (auto k = toSchedule.begin();k != toSchedule.end();k++) {
				 totalCompleted.insert(*k);
			 }

			 //removing redundant gates
			 for (auto k = currentGates.begin();k != currentGates.end();) {
				 if (totalCompleted.find(*k) != totalCompleted.end()) {
					 k = currentGates.erase(k);
				 }
				 else {
					 k++;
				 }
			 }
			 
			 //we do scheduling for the initial maximalMapper
			 /*
			 while (forSchedule.size() > 0) {
				 toSchedule.insert(forSchedule.front());
				 forSchedule.pop();
			 }
			 toSchedule = get<3>((*i));*/

			 //doing scheduling just like in schedule method for optimal
			 int numToMap = toSchedule.size();
			 while (finished.size() < numToMap) {
				 for (auto i = toSchedule.begin();i != toSchedule.end();) {

					 GateNode* toConsider = *i;
					 if ((*i)->control != -1) {
						 //2 qubit gate
						 //if the parents are already in finish set, or neither in 	
						 if ((toConsider->targetParent == NULL ||
							 finished.find(toConsider->targetParent) != finished.end() ||
							 (toSchedule.find(toConsider->targetParent) == toSchedule.end() && inProgress.find(toConsider->targetParent) == inProgress.end())) &&
							 ((toConsider->controlParent == NULL) || finished.find(toConsider->controlParent) != finished.end() ||
								 (toSchedule.find(toConsider->controlParent) == toSchedule.end() && inProgress.find(toConsider->controlParent) == inProgress.end()))) {
							 //then both parents are satisfied, lets try to map these qubits
							 if (logicalSchedule[toConsider->target].second == -1 && logicalSchedule[toConsider->control].second == -1) {
								 //then we can schedule this gate for both qubits
								 //printf("LATENCYDOUBLE %d\n", latencyDouble);
								 logicalSchedule[toConsider->target] = make_pair(toConsider, latencyDouble);
								 logicalSchedule[toConsider->control] = make_pair(toConsider, latencyDouble);
								 //removing this gate from schedule list
								 inProgress.insert(toConsider);
								 //removing gate from toSchedule
								 i = toSchedule.erase(i);
								 //incrementing gate number
								 numGates++;

								 //debug--> making sure that this is valid
								 int firstPhys = toTransform[toConsider->control];
								 int secondPhys = toTransform[toConsider->target];
								 if (architectureEdges[firstPhys].find(secondPhys) == architectureEdges[firstPhys].end()) {

									 cout << " SOMETHING WENT TERRIBLY WRONG IN SCHEDULING A QUBIT CNOT!\n";
								 }

							 }
							 else {
								 //moving on 
								 i++;
							 }
						 }
						 else {
							 //move on
							 i++;
						 }
					 }
					 else {
						 //1 qubit gate	
						 //check first if parent is being used or not
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
								 i = toSchedule.erase(i);

								 //incrementing the gate number
								 numGates++;
							 }
							 else {
								 i++;
							 }
						 }
						 else {
							 //move on
							 i++;
						 }
					 }
				 }
				 //simulating a cycle here
				for (int z = 0; z < numLogical;z++) {
					if (logicalSchedule[z].second > 0) {
						//decrementing
						GateNode* work = logicalSchedule[z].first;
						int cyclesLeft = logicalSchedule[z].second;
						cyclesLeft--;
						 if (cyclesLeft == 0) {
							 //then this gate is done executing, we should remove it from the inProgress set and from logicalSchedule
							 inProgress.erase(work);
							 finished.insert(work);
							 //if this gate is a swap gate, then we transform the mapping
							 //need to make both logical and control qubits the dummy
							 logicalSchedule[work->target] = make_pair(dummy, -1);
							 if (work->control != -1) {
								 logicalSchedule[work->control] = make_pair(dummy, -1);
							 }
							 if (work->name == "swp") {
								 numSwaps++;
								 numGates++;
								 //we do the swap change
								 //debug --> checking if swap is valid or not
								 int physOne = toTransform[work->control];
								 int physTwo = toTransform[work->target];
								 
								 if (architectureEdges[physOne].find(physTwo) == architectureEdges[physOne].end()) {
									cout << "SOMETHING WENT TERRIBLY WRONG WITH MAPPING FOR SWAP!\n";
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
						 else {
							 //update the cycle value 
							 logicalSchedule[z] = make_pair(work, cyclesLeft);
						 }
						
					}
				}
				//this simulates a single cycle
				numCycles++;
	
			 }

			 //we are good here, now we empty the sets
			 toSchedule.clear();
			 inProgress.clear();
			 finished.clear();


			 continue;
		 }

		 if (targetMapping.size() == 0) {
			 //generating closest mapping to the last mapping
			 //first cropping the architecture to our initial mapping and then updating distances
			 architectureEdges = getCroppedArchitecture(initialMapping, architectureEdges);

			 distances = getDistances(architectureEdges);
			 
			 vector<tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>>> res = maximalMapper(architectureEdges, currentGates, true,totalCompleted);
			 //picking the "closest" mapping to our last mapping (transformed Mapping)
			 int cost = -1;
			 map<int, int> closest;
			 set<GateNode*> endGates;
			 for (auto k = res.begin();k != res.end();k++) {
				 int stitchingCost = partialMappingStitchingCost(toTransform, get<0>(*k), distances, architectureEdges);
				 if (cost == -1 || stitchingCost < cost) {
					 cost = stitchingCost;
					 closest = get<0>(*k);
					 endGates = get<2>(*k);
					 toSchedule = get<3>(*k);
				}
			 }

			 //adding scheduled gates to totalCompletedSet
			 for (auto k = toSchedule.begin();k != toSchedule.end();k++) {
				 totalCompleted.insert(*k);
			 }


			currentGates = endGates;
			for (auto k = currentGates.begin();k != currentGates.end();) {
				 if (totalCompleted.find(*k) != totalCompleted.end()) {
					 k = currentGates.erase(k);
				 }
				 else {
					 k++;
				 }
			 }

			 //getting swaps needed to stitch
			 pair<map<int,int>,vector<pair<int,int>>> swaps = stitchMappings(make_pair(toTransform, set<int>()), make_pair(closest, set<int>()), architectureEdges, distances);
			 
			 //parallelizing swaps 
			 vector<set<pair<int, int>>> parallelizedSwaps = parallelizeSwapPathway(swaps.second);
			 //tacking on swaps to end of partition basically, and we will print and apply the swaps 1 by 1
			 numCycles += parallelizedSwaps.size() * latencySwap;
			 for (auto i = parallelizedSwaps.begin();i != parallelizedSwaps.end();i++) {
				 for (auto j = (*i).begin();j != (*i).end();j++) {
					 //printing the swap, then applying it 
					 numSwaps++;
					 numGates++;
					 cout << "swp q[" << toTransform[j->first] << "], q[" << toTransform[j->second] << "] ; // between logical q[" << j->first << "], q[" << j->second << "]\n";
					 //performing the swap
					 int tmp = toTransform[j->first];
					 toTransform[j->first] = toTransform[j->second];
					 toTransform[j->second] = tmp;
				 }
			 }

			 //we are done applying swaps
			 //now we get the next queue to print
			 //scheduling these gates
			 /*
			 queue<GateNode*> satisfiedGates = getSatisfiedGates(currentGates, endGates);
			 currentGates = endGates;
			 while (satisfiedGates.size() > 0) {
				//do the same thing as above
				 toSchedule.insert(satisfiedGates.front());
				 satisfiedGates.pop();
			 }*/

			 int numToMap = toSchedule.size();
			 while (finished.size() < numToMap) {
				 for (auto i = toSchedule.begin();i != toSchedule.end();) {

					 GateNode* toConsider = *i;
					 if ((*i)->control != -1) {
						 //2 qubit gate
						 //if the parents are already in finish set, or neither in 	
						 if ((toConsider->targetParent == NULL ||
							 finished.find(toConsider->targetParent) != finished.end() ||
							 (toSchedule.find(toConsider->targetParent) == toSchedule.end() && inProgress.find(toConsider->targetParent) == inProgress.end())) &&
							 ((toConsider->controlParent == NULL) || finished.find(toConsider->controlParent) != finished.end() ||
								 (toSchedule.find(toConsider->controlParent) == toSchedule.end() && inProgress.find(toConsider->controlParent) == inProgress.end()))) {
							 //then both parents are satisfied, lets try to map these qubits
							 if (logicalSchedule[toConsider->target].second == -1 && logicalSchedule[toConsider->control].second == -1) {
								 //then we can schedule this gate for both qubits
								 logicalSchedule[toConsider->target] = make_pair(toConsider, latencyDouble);
								 logicalSchedule[toConsider->control] = make_pair(toConsider, latencyDouble);
								 //removing this gate from schedule list
								 inProgress.insert(toConsider);
								 //removing gate from toSchedule
								 i = toSchedule.erase(i);
								 //incrementing gate number
								 numGates++;

								 //debug--> making sure that this is valid
								 int firstPhys = toTransform[toConsider->control];
								 int secondPhys = toTransform[toConsider->target];
								 if (architectureEdges[firstPhys].find(secondPhys) == architectureEdges[firstPhys].end()) {

									 cout << " SOMETHING WENT TERRIBLY WRONG IN SCHEDULING A QUBIT CNOT!\n";
								 }

							 }
							 else {
								 //moving on 
								 i++;
							 }
						 }
						 else {
							 //move on 
							 i++;
						 }
					 }
					 else {
						 //1 qubit gate	
						 //check first if parent is being used or not
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
								 i = toSchedule.erase(i);

								 //incrementing the gate number
								 numGates++;
							 }
							 else {
								 i++;
							 }
						 }
						 else {
							 //move on
							 i++;
						 }
					 }
				 }
				 //simulating a cycle here
				for (int z = 0; z < numLogical;z++) {
					if (logicalSchedule[z].second > 0) {
						//decrementing
						GateNode* work = logicalSchedule[z].first;
						int cyclesLeft = logicalSchedule[z].second;
						cyclesLeft--;
						 if (cyclesLeft == 0) {
							 //then this gate is done executing, we should remove it from the inProgress set and from logicalSchedule
							 inProgress.erase(work);
							 finished.insert(work);
							 //if this gate is a swap gate, then we transform the mapping
							 //need to make both logical and control qubits the dummy
							 logicalSchedule[work->target] = make_pair(dummy, -1);
							 if (work->control != -1) {
								 logicalSchedule[work->control] = make_pair(dummy, -1);
							 }
							 if (work->name == "swp") {
								 numSwaps++;
								 //we do the swap change
								 //debug --> checking if swap is valid or not
								 int physOne = toTransform[work->control];
								 int physTwo = toTransform[work->target];
								 
								 if (architectureEdges[physOne].find(physTwo) == architectureEdges[physOne].end()) {
									cout << "SOMETHING WENT TERRIBLY WRONG WITH MAPPING FOR SWAP!\n";
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
						 else {
							 logicalSchedule[z] = make_pair(work, cyclesLeft);
						 }
						
					}
				}
				//this simulates a single cycle
				numCycles++;

			 }


			 //setting second mapping back to nothing (since now transformed mapping handles it)
			 targetMapping = map<int, int>();
			 //clearing the sets
			 toSchedule.clear();
			 inProgress.clear();
			 finished.clear();
		 }
	 }
	 cout << " NUMBER GATES: " << numGates <<"\n";
	 cout << " NUMBER SWAPS: " << numSwaps << "\n";
	 cout << " NUMBER CYCLES: " << numCycles << "\n";

	 //printing out initial mapping
	 cout << "INITIAL MAPPING!\n";
	 for (auto i = initialMapping.begin();i != initialMapping.end();i++) {
		 cout << "Logical q[" << i->first << "] mapped to physical Q[" << i->second<<"] \n";
	 }
}

 //idea here is that as we go, we can print so that we do not take up space
	//this is my fast mapper, so we do not have to keep all the memory required until we get an optimal solution
 void parallelSchedulingFastCircuitBuilder(map<int,set<int>> originalArchitectureEdges, set<GateNode*> startGates) {
	 //will eventually need to print this out
	 int numSwaps = 0;
	 int numGates = 0;
	 int numCycles = 0;
	 map<pair<int, int>, int> distances = getDistances(originalArchitectureEdges);
		
	 queue<GateNode*> forSchedule;
	 set<GateNode*> toSchedule;
	 set<GateNode*> inProgress;
	 set<GateNode*> finished;
	 vector<pair<GateNode*, int>> logicalSchedule;
	 GateNode* dummy = new GateNode();
	 for (int i = 0;i < numLogical;i++) {
		 //-1 means qubit is free to participate in some mapping
		 logicalSchedule.push_back( make_pair(dummy, -1));
	 }
	 //maintaining our initial mapping (every qubit may not be mapped at the start)
	 map<int, int> initialMapping;
	 //keeping track of our current state
	 map<int, int> toTransform;
	 //keeping track of the target mapping for swaps
	 map<int, int> targetMapping;
	 set<GateNode*> currentGates = startGates;
	 map<int, set<int>> architectureEdges = originalArchitectureEdges;

	 while (currentGates.size() > 0) {

		 if (initialMapping.size() == 0) {
			 //need to get initial maximal mapping
			 vector<tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>>> partialMappings = maximalMapper(architectureEdges, currentGates, true,set<GateNode*>());
			 //picking lowest cost from these
			 map<int, int> bestMapping;
			 set<GateNode*> endGates;
			 int cost = -1;
			 for (auto i = partialMappings.begin();i != partialMappings.end();i++) {
				 int mapCost = partialMappingCost(get<0>(*i), get<2>(*i), distances, architectureEdges);
				 if (cost == -1 || mapCost < cost) {
					 bestMapping = get<0>(*i);
					 cost = mapCost;
					 endGates = get<2>(*i);
				 }
			 }

			 //filling the initial mapping
			 map<int, int> fullInitialMapping = pickBestFilledMapping(bestMapping, architectureEdges, endGates, distances);
			 initialMapping = fullInitialMapping;
			 
			 toTransform = fullInitialMapping;
			 forSchedule = getSatisfiedGates(currentGates, endGates);
			 currentGates = endGates;
			 
			 //we do scheduling for the initial maximalMapper
			 while (forSchedule.size() > 0) {
				 toSchedule.insert(forSchedule.front());
				 forSchedule.pop();
			 }

			 //doing scheduling just like in schedule method for optimal
			 int numToMap = toSchedule.size();
			 while (finished.size() < numToMap) {
				 for (auto i = toSchedule.begin();i != toSchedule.end();) {

					 GateNode* toConsider = *i;
					 if ((*i)->control != -1) {
						 //2 qubit gate
						 //if the parents are already in finish set, or neither in 	
						 if ((toConsider->targetParent == NULL ||
							 finished.find(toConsider->targetParent) != finished.end() ||
							 (toSchedule.find(toConsider->targetParent) == toSchedule.end() && inProgress.find(toConsider->targetParent) == inProgress.end())) &&
							 ((toConsider->controlParent == NULL) || finished.find(toConsider->controlParent) != finished.end() ||
								 (toSchedule.find(toConsider->controlParent) == toSchedule.end() && inProgress.find(toConsider->controlParent) == inProgress.end()))) {
							 //then both parents are satisfied, lets try to map these qubits
							 if (logicalSchedule[toConsider->target].second == -1 && logicalSchedule[toConsider->control].second == -1) {
								 //then we can schedule this gate for both qubits
								 //printf("LATENCYDOUBLE %d\n", latencyDouble);
								 logicalSchedule[toConsider->target] = make_pair(toConsider, latencyDouble);
								 logicalSchedule[toConsider->control] = make_pair(toConsider, latencyDouble);
								 //removing this gate from schedule list
								 inProgress.insert(toConsider);
								 //removing gate from toSchedule
								 i = toSchedule.erase(i);
								 //incrementing gate number
								 numGates++;

								 //debug--> making sure that this is valid
								 int firstPhys = toTransform[toConsider->control];
								 int secondPhys = toTransform[toConsider->target];
								 if (architectureEdges[firstPhys].find(secondPhys) == architectureEdges[firstPhys].end()) {

									 cout << " SOMETHING WENT TERRIBLY WRONG IN SCHEDULING A QUBIT CNOT!\n";
								 }

							 }
							 else {
								 //moving on 
								 i++;
							 }
						 }
						 else {
							 //move on
							 i++;
						 }
					 }
					 else {
						 //1 qubit gate	
						 //check first if parent is being used or not
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
								 i = toSchedule.erase(i);

								 //incrementing the gate number
								 numGates++;
							 }
							 else {
								 i++;
							 }
						 }
						 else {
							 //move on
							 i++;
						 }
					 }
				 }
				 //simulating a cycle here
				for (int z = 0; z < numLogical;z++) {
					if (logicalSchedule[z].second > 0) {
						//decrementing
						GateNode* work = logicalSchedule[z].first;
						int cyclesLeft = logicalSchedule[z].second;
						cyclesLeft--;
						 if (cyclesLeft == 0) {
							 //then this gate is done executing, we should remove it from the inProgress set and from logicalSchedule
							 inProgress.erase(work);
							 finished.insert(work);
							 //if this gate is a swap gate, then we transform the mapping
							 //need to make both logical and control qubits the dummy
							 logicalSchedule[work->target] = make_pair(dummy, -1);
							 if (work->control != -1) {
								 logicalSchedule[work->control] = make_pair(dummy, -1);
							 }
							 if (work->name == "swp") {
								 numSwaps++;
								 numGates++;
								 //we do the swap change
								 //debug --> checking if swap is valid or not
								 int physOne = toTransform[work->control];
								 int physTwo = toTransform[work->target];
								 
								 if (architectureEdges[physOne].find(physTwo) == architectureEdges[physOne].end()) {
									cout << "SOMETHING WENT TERRIBLY WRONG WITH MAPPING FOR SWAP!\n";
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
						 else {
							 //update the cycle value 
							 logicalSchedule[z] = make_pair(work, cyclesLeft);
						 }
						
					}
				}
				//this simulates a single cycle
				numCycles++;
	
			 }

			 //we are good here, now we empty the sets
			 toSchedule.clear();
			 inProgress.clear();
			 finished.clear();


			 continue;
		 }

		 if (targetMapping.size() == 0) {
			 //generating closest mapping to the last mapping
			 //first cropping the architecture to our initial mapping and then updating distances
			 architectureEdges = getCroppedArchitecture(initialMapping, architectureEdges);

			 distances = getDistances(architectureEdges);
			 
			 vector<tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>>> res = maximalMapper(architectureEdges, currentGates, true,set<GateNode*>());
			 //picking the "closest" mapping to our last mapping (transformed Mapping)
			 int cost = -1;
			 map<int, int> closest;
			 set<GateNode*> endGates;
			 for (auto k = res.begin();k != res.end();k++) {
				 int stitchingCost = partialMappingStitchingCost(toTransform, get<0>(*k), distances, architectureEdges);
				 if (cost == -1 || stitchingCost < cost) {
					 cost = stitchingCost;
					 closest = get<0>(*k);
					 endGates = get<2>(*k);
				}
			 }
			

			 //getting swaps needed to stitch
			 pair<map<int,int>,vector<pair<int,int>>> swaps = stitchMappings(make_pair(toTransform, set<int>()), make_pair(closest, set<int>()), architectureEdges, distances);
			 
			 //parallelizing swaps 
			 vector<set<pair<int, int>>> parallelizedSwaps = parallelizeSwapPathway(swaps.second);
			 //tacking on swaps to end of partition basically, and we will print and apply the swaps 1 by 1
			 numCycles += parallelizedSwaps.size() * latencySwap;
			 for (auto i = parallelizedSwaps.begin();i != parallelizedSwaps.end();i++) {
				 for (auto j = (*i).begin();j != (*i).end();j++) {
					 //printing the swap, then applying it 
					 numSwaps++;
					 numGates++;
					 cout << "swp q[" << toTransform[j->first] << "], q[" << toTransform[j->second] << "] ; // between logical q[" << j->first << "], q[" << j->second << "]\n";
					 //performing the swap
					 int tmp = toTransform[j->first];
					 toTransform[j->first] = toTransform[j->second];
					 toTransform[j->second] = tmp;
				 }
			 }

			 //we are done applying swaps
			 //now we get the next queue to print
			 //scheduling these gates
			 queue<GateNode*> satisfiedGates = getSatisfiedGates(currentGates, endGates);
			 currentGates = endGates;
			 while (satisfiedGates.size() > 0) {
				//do the same thing as above
				 toSchedule.insert(satisfiedGates.front());
				 satisfiedGates.pop();
			 }

			 int numToMap = toSchedule.size();
			 while (finished.size() < numToMap) {
				 for (auto i = toSchedule.begin();i != toSchedule.end();) {

					 GateNode* toConsider = *i;
					 if ((*i)->control != -1) {
						 //2 qubit gate
						 //if the parents are already in finish set, or neither in 	
						 if ((toConsider->targetParent == NULL ||
							 finished.find(toConsider->targetParent) != finished.end() ||
							 (toSchedule.find(toConsider->targetParent) == toSchedule.end() && inProgress.find(toConsider->targetParent) == inProgress.end())) &&
							 ((toConsider->controlParent == NULL) || finished.find(toConsider->controlParent) != finished.end() ||
								 (toSchedule.find(toConsider->controlParent) == toSchedule.end() && inProgress.find(toConsider->controlParent) == inProgress.end()))) {
							 //then both parents are satisfied, lets try to map these qubits
							 if (logicalSchedule[toConsider->target].second == -1 && logicalSchedule[toConsider->control].second == -1) {
								 //then we can schedule this gate for both qubits
								 logicalSchedule[toConsider->target] = make_pair(toConsider, latencyDouble);
								 logicalSchedule[toConsider->control] = make_pair(toConsider, latencyDouble);
								 //removing this gate from schedule list
								 inProgress.insert(toConsider);
								 //removing gate from toSchedule
								 i = toSchedule.erase(i);
								 //incrementing gate number
								 numGates++;

								 //debug--> making sure that this is valid
								 int firstPhys = toTransform[toConsider->control];
								 int secondPhys = toTransform[toConsider->target];
								 if (architectureEdges[firstPhys].find(secondPhys) == architectureEdges[firstPhys].end()) {

									 cout << " SOMETHING WENT TERRIBLY WRONG IN SCHEDULING A QUBIT CNOT!\n";
								 }

							 }
							 else {
								 //moving on 
								 i++;
							 }
						 }
						 else {
							 //move on 
							 i++;
						 }
					 }
					 else {
						 //1 qubit gate	
						 //check first if parent is being used or not
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
								 i = toSchedule.erase(i);

								 //incrementing the gate number
								 numGates++;
							 }
							 else {
								 i++;
							 }
						 }
						 else {
							 //move on
							 i++;
						 }
					 }
				 }
				 //simulating a cycle here
				for (int z = 0; z < numLogical;z++) {
					if (logicalSchedule[z].second > 0) {
						//decrementing
						GateNode* work = logicalSchedule[z].first;
						int cyclesLeft = logicalSchedule[z].second;
						cyclesLeft--;
						 if (cyclesLeft == 0) {
							 //then this gate is done executing, we should remove it from the inProgress set and from logicalSchedule
							 inProgress.erase(work);
							 finished.insert(work);
							 //if this gate is a swap gate, then we transform the mapping
							 //need to make both logical and control qubits the dummy
							 logicalSchedule[work->target] = make_pair(dummy, -1);
							 if (work->control != -1) {
								 logicalSchedule[work->control] = make_pair(dummy, -1);
							 }
							 if (work->name == "swp") {
								 numSwaps++;
								 //we do the swap change
								 //debug --> checking if swap is valid or not
								 int physOne = toTransform[work->control];
								 int physTwo = toTransform[work->target];
								 
								 if (architectureEdges[physOne].find(physTwo) == architectureEdges[physOne].end()) {
									cout << "SOMETHING WENT TERRIBLY WRONG WITH MAPPING FOR SWAP!\n";
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
						 else {
							 logicalSchedule[z] = make_pair(work, cyclesLeft);
						 }
						
					}
				}
				//this simulates a single cycle
				numCycles++;

			 }


			 //setting second mapping back to nothing (since now transformed mapping handles it)
			 targetMapping = map<int, int>();
			 //clearing the sets
			 toSchedule.clear();
			 inProgress.clear();
			 finished.clear();
		 }
	 }
	 cout << " NUMBER GATES: " << numGates <<"\n";
	 cout << " NUMBER SWAPS: " << numSwaps << "\n";
	 cout << " NUMBER CYCLES: " << numCycles << "\n";

	 //printing out initial mapping
	 cout << "INITIAL MAPPING!\n";
	 for (auto i = initialMapping.begin();i != initialMapping.end();i++) {
		 cout << "Logical q[" << i->first << "] mapped to physical Q[" << i->second<<"] \n";
	 }
	
 }
/*
 //my faster version of the optimal solver
 void betterLazySwapCircuitBuilder(map<int,set<int>> originalArchitectureEdges, set<GateNode*> startGates) {
 //since our initial mapping varies until everything is mapped, we have to print out gates at the end
 //but the process is this: find the best initial maximal mapping (maximal Mapping cost)
	 //gates satisfied for every step of the maximal mapping
	 //distance matrix to be updated at each step
	 map<pair<int, int>, int> distances = getDistances(originalArchitectureEdges);
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
		 vector<tuple<map<int,int>,set<int>,set<GateNode*>,set<GateNode*>>> partialMapping=maximalMapper(architectureEdges, currentGates,true,set<GateNode*>());

		 
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

			 //updating the last mapping
			 set<int> toContinue = lastMapping.second;
			 lastMapping = make_pair(filledMapping, toContinue);

			 //setting our initial mapping to last mapping
			 initialMapping = lastMapping.first;

			 //updating architecture edges
			 architectureEdges = getCroppedArchitecture(filledMapping, architectureEdges);

			 //updating distance matrix for our new architecture
			 distances = getDistances(architectureEdges);

			 //adding gates to the print queue
			 queue<GateNode*> toPrint = getSatisfiedGates(currentGates, endGates);

			 //updating currentGates
			 currentGates = endGates;

			//pushing print queue to vector
			maximalMappingSatisfiedGates.push_back(toPrint);
			continue;
		 }
		 else {
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
			 //now we find the list of swaps to take
			 pair<map<int, int>, vector<pair<int, int>>> swaps = stitchMappings(lastMapping, bestMapping, architectureEdges,distances);

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
		 }


	 }
	 
	 //we call schedule here
	
	 scheduleCircuit(maximalMappingSatisfiedGates,levelSwaps,initialMapping,architectureEdges);

	 //printing the initial mapping
	 cout << "INITIAL MAPPING!: \n";
	 for (auto z = initialMapping.begin();z != initialMapping.end();z++) {
		 cout << "Logical qubit q" << z->first << " mapped to physical qubit Q" << z->second <<"\n";
	 }
 }*/



/*
 //function to test my maximal mapper
	//it should return mappings for every maximal portion of the circuit that can be executed without swaps
 queue<map<int, int>> getMaximalMappings(set<GateNode*> startGates,map<int,set<int>> architectureEdges) {
	 set<GateNode*> nextGates = startGates;

	 queue<map<int, int>> resultQueue;
	 int count = 1;
	 do {
		vector <tuple<map<int, int>, set<int>, set<GateNode*>,set<GateNode*>> > result = maximalMapper(architectureEdges,nextGates,false,set<GateNode*>());
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
}*/

 
