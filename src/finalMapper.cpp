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
#include <algorithm>
#include <unordered_map>
#include <list>
#include <queue>
#include <map>
#include <deque>

//if smart is defined, we will use a distance matrix instead of going through bfs every time
int SMART = 1;
//if latency is 1, we will try and gracefully insert swaps between maximal partitions
int LATENCY= 1;

vector <tuple<map<int, int>, set<int>, set<GateNode*>> > maximalMapper(map<int, int> architectureNodeDegrees, multimap<int, int> architectureEdges, set<GateNode*> startSet);
pair<int, int> addGateEdge(GateNode* chosen, multimap<int, int>circuitEdges, map<int, int>circuitNodeDegrees);
queue<pair<int, int>> pruneArchitectureEdges(multimap<int, int>architectureEdges, map<int, int>architectureNodeDegrees, int circuitMinDegree);
bool wentTooFar(map<int, int>circuitNodeDegrees, map<int, int>architectureNodeDegrees, int circuitMaxDegree, int architectureMaxDegree);
void rollBack(queue<pair<int, int>>changes, multimap<int, int>architectureEdges, map<int, int>architectureNodeDegrees);


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
	//debug
	cout << "PRINTING COUPLING\n";
	for (auto i = couplings.begin();i != couplings.end();i++) {
		printf("%d connected to %d\n", (*i).first, (*i).second);
	}

	//creating a multimap for quick lookup of architecture edges
	multimap<int, int> architectureEdges ;
	map<int, int> architectureNodeDegrees;

	//building edges and degrees
	for (auto i = couplings.begin();i != couplings.end();i++) {
		architectureEdges.insert(make_pair(i->first, i->second));
		architectureEdges.insert(make_pair(i->second, i->first));

		//incrementing degrees
		if (architectureNodeDegrees.find(i->first) == architectureNodeDegrees.end()) {
			architectureNodeDegrees[i->first] = 1;
		}
		else {
			architectureNodeDegrees[i->first]++;
		}

		if (architectureNodeDegrees.find(i->second) == architectureNodeDegrees.end()) {
			architectureNodeDegrees[i->second] = 1;
		}
		else {
			architectureNodeDegrees[i->second]++;
		}
	}

	cout << "Printing architecture edges!\n";
	for (auto i = architectureEdges.begin();i != architectureEdges.end();i++) {
		printf("Physical Edge between physical qubit %d and physical qubit %d\n", i->first, i->second);
	}

	cout << "Printing architecture degrees!\n";
	for (auto i = architectureNodeDegrees.begin();i != architectureNodeDegrees.end();i++) {
		printf("Physical qubit %d has degree %d\n", i->first, i->second);
	}

			
	//right now calling maximal mapper once for perfect mapping problem
	vector <tuple<map<int, int>, set<int>, set<GateNode*>>> result = maximalMapper(architectureNodeDegrees, architectureEdges, firstGates);
	tuple<map<int, int>, set<int>, set<GateNode*>> anyMapping = result[0];

	if (get<2>(anyMapping).size() == 0) {
		//then we found a perfect mapping
		map<int, int> perfectMapping = get<0>(anyMapping);
		cout << "Found a Perfect Initial Mapping!\n";
		printf("Location of qubits: ");
		for (int i = 0;i < numLogicalQubits;i++) {
			printf("%d, ", perfectMapping.find(i)->second);
		}
		printf("\n");
	}
	else {
		//then we did not find a perfect mapping
		cout << "No Perfect Initial Mapping Found!\n";
	}

	//Exit the program:
	return 0;
}

//maximal mapper that returns matchings for the next maximal portion of the dependence graph
	//returns a pair
		//the first element is a list of all maximal matchings
		//the second element is jsut a set of mapped architecture nodes --> for quick lookup
		//the third element is the queue of gatenodes which could not be matched and should be considered the "start" for the next iteration
		//this is returned so we can easily start the next iteration of mapping
vector <tuple<map<int, int>, set<int>, set<GateNode*>> > maximalMapper(map<int, int> architectureNodeDegrees, multimap<int, int> architectureEdges, set<GateNode*> startSet) {
	//we first grab a copy of the coupling graph and associated nodeDegree vector
	multimap<int, int> copyArchitectureEdges = architectureEdges;
	map<int,int> copyArchitectureNodeDegrees = architectureNodeDegrees;

	//allocating graph structure to add level by level in the dependence graph
	multimap<int, int> circuitEdges ;
	map<int,int> circuitNodeDegrees ;


	//keeping track of considered gates (we do not want to consider gates twice)
	set<GateNode*> closedSet;
	//this is a queue to insert accounted gates into, will come into play when we actually start mapping
	deque<GateNode*> mappingQueue;
	//vector<GateNode*> mappingQueue;
	//set<GateNode*> mappingSet;
	set<GateNode*> couldntMap;

		
	queue<GateNode*> startGates;
	set<GateNode*> start = startSet;
	while (start.size() > 0) {
		for (auto i = start.begin();i != start.end();i++) {
			
			if (start.find((*i)->controlParent) == start.end() && start.find((*i)->targetParent) == start.end()) {
				//then we add it into queue to get a level order
				startGates.push((*i));
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
		pair<int,int> minMaxDegrees = addGateEdge(chosen, circuitEdges,circuitNodeDegrees);
		queue<pair<int,int>> changes = pruneArchitectureEdges(copyArchitectureEdges,copyArchitectureNodeDegrees,minMaxDegrees.first);
		
		
		int circuitMaxDegree = minMaxDegrees.second;
		//suppose we went too far
			//then we rollback and continue with the queue
			// keep doing this until queue is emtpy --> that is our maximal partition
		
		//getting max degree of the architecture nodes
		int architectureMaxDegree = 0;
		for (auto i = copyArchitectureNodeDegrees.begin();i != copyArchitectureNodeDegrees.end();i++) {
			if (i->second > architectureMaxDegree) {
				architectureMaxDegree = i->second;
			}
		}
		if (wentTooFar(circuitNodeDegrees,copyArchitectureNodeDegrees,circuitMaxDegree,architectureMaxDegree)) {
			//rollback 
			rollBack(changes, architectureEdges, architectureNodeDegrees);
			//undoing circuit node add and edges/degree

			//we couldnt map the given gate, so we will push it to a list of nodes 
			couldntMap.insert(chosen);
			//we do not want to consider this gate again
			closedSet.insert(chosen);
		}
		else {
			//then we can update (add children of the block to the queue
			if (chosen->controlChild != NULL ) {
				//then we need to check if target parent is in the closedSet or not for insertion (so we enforce level order)
				if (closedSet.find((chosen->controlChild)->targetParent) != closedSet.end()) {
					//then we can add this child to the queue
					startGates.push(chosen->controlChild);
				}
			}

			if (chosen->targetChild != NULL ) {
				//then we need to check if control parent is in the closedSet or not for insertion (so we enforce level order)
				if (closedSet.find((chosen->targetChild)->controlParent) != closedSet.end()) {
					//then we can add this child to the queue
					startGates.push(chosen->targetChild);
				}
			}

			//adding the chosen gate to the mapping queue
			mappingQueue.push_back(chosen);
			closedSet.insert(chosen);
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
		
	//our stack is of the form <(control, pairing),(target,pairing)>
	//stack < pair<pair<int, int>,pair<int,int>>> maximalMappingStack;
	//using a queue to go gate by gate and extend our maximal mappings
	queue < tuple<map<int, int>,set<int>, set<GateNode*>>> maximalMappingQueue;
	//keeping track of where we are in the mapping
	set<GateNode*> otherClosedSet; 
	while (mappingQueue.size() > 0) {
		GateNode* toSatisfy = mappingQueue.front();
		mappingQueue.pop_front();
		//preserving order for return values
		if (maximalMappingQueue.empty()) {
			//we do the first iteration
			GateNode* firstGate = toSatisfy;
			otherClosedSet.insert(firstGate);
			int control = firstGate->control;
			int target = firstGate->target;
			auto iter = copyArchitectureNodeDegrees.begin();
			while (iter != copyArchitectureNodeDegrees.end()) {
				if (iter->second >= circuitNodeDegrees[control]) {
					//then we found a physical qubit that can accomodate the logical qubit "control" in this maximal mapping
					//we need to check the neighbors of this qubit for possible "target" mappings
					auto edgeRange = copyArchitectureEdges.equal_range(iter->first);
					auto edgeIter = edgeRange.first;
					while (edgeIter != edgeRange.second) {
						if (copyArchitectureNodeDegrees[(edgeIter->second)] >= circuitNodeDegrees[target]) {
							//then this is a possible mapping for the target logical qubit
							map<int, int> suggestedInitialMapping;
							set<int> mappedArchitectureQubits;
							suggestedInitialMapping.insert(make_pair(control, iter->first));
							mappedArchitectureQubits.insert(iter->first);
							suggestedInitialMapping.insert(make_pair(target, edgeIter->second));
							mappedArchitectureQubits.insert(edgeIter->second);
							set<GateNode*> currentlySatisfied;
							currentlySatisfied.insert(firstGate);
							//pushing the mapping and satisfied gate to the queue 	
							maximalMappingQueue.push(make_tuple(suggestedInitialMapping,mappedArchitectureQubits, currentlySatisfied));
						}
						edgeIter++;
					}
				}
				iter++;
			}
			continue;
		} 	

		//we go through all mappings in the queue and see if we can propogate them to the next gate
			//we can easily check if parent gates are satisfied in the current mapping or not with the extra set given
		queue < tuple<map<int, int>, set<int>,set<GateNode*>>> resultMaximalMappingQueue;
		while (maximalMappingQueue.size() > 0) {
			tuple<map<int, int>, set<int>, set<GateNode*>> queueEntry = maximalMappingQueue.front();
			maximalMappingQueue.pop();
			//seeing if we can propagate the mapping to the next gate, and inserting it in the queue
			map<int, int> currMapping = get<0>(queueEntry);
			set<int> mappedArchQubits = get<1>(queueEntry);
			set<GateNode*> satisfiedGates = get<2>(queueEntry);
			
			//firstly checking if the parent gates part of the prior mapping sequences and if so, they should be mapped in the parents
				//if not, then we cannot map this gate
			if ((otherClosedSet.find(toSatisfy->controlParent) != otherClosedSet.end() && satisfiedGates.find(toSatisfy->controlParent) == satisfiedGates.end()
				) || (otherClosedSet.find(toSatisfy->targetParent) != otherClosedSet.end() && satisfiedGates.find(toSatisfy->targetParent) == satisfiedGates.end())) {
				//then the parent gates are not satisfied, we cannot propogate (extend) this mapping to the next gate
				resultMaximalMappingQueue.push(queueEntry);
				continue;
			}

			//dependencies are resolved if we reached here
				//we need to check 3 cases
				//1) both qubits are mapped
				//2) only 1 of the two qubits are mapped
				//3) none of the qubits are mapped
			if (currMapping.find(toSatisfy->control) != currMapping.end() && currMapping.find(toSatisfy->target) != currMapping.end()) {
				//then both are already mapped, we just check if they satisfy the new gate to satisfy
				//namely the architecture nodes they are mapped to should have an edge between them
				auto range = copyArchitectureEdges.equal_range(currMapping[toSatisfy->control]);
				bool found = false;
				for (auto i = range.first;i != range.second;i++) {
					if (i->second == currMapping[toSatisfy->target]) {
						//then we can propogate this mapping to next gate level 
						found = true;
						break;
					}
				}
				if (found) {
					satisfiedGates.insert(toSatisfy);
					resultMaximalMappingQueue.push(make_tuple(currMapping,mappedArchQubits, satisfiedGates));
				}
				else {
					//no change in satisfied gates or the mapping
					resultMaximalMappingQueue.push(queueEntry);
				}
			}
			else if (currMapping.find(toSatisfy->control) != currMapping.end() || currMapping.find(toSatisfy->target) != currMapping.end()) {
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
				auto range = copyArchitectureEdges.equal_range(currMapping[mappedQubit]);
				for (auto i = range.first;i != range.second;i++) {
					if (copyArchitectureNodeDegrees[i->second] >= circuitNodeDegrees[unmappedQubit] && mappedArchQubits.find(i->second)==mappedArchQubits.end()) {
						//then this is a possible mapping to propogate
						didPropagate = true;
						map<int, int> newMapping = currMapping;
						set<int> archMapping = mappedArchQubits;
						archMapping.insert(i->second);
						newMapping.insert(make_pair(unmappedQubit, i->second));
						set<GateNode*> newSatisfied = satisfiedGates;
						newSatisfied.insert(toSatisfy);
						resultMaximalMappingQueue.push(make_tuple(newMapping,archMapping, newSatisfied));
					}
				}
				if (!didPropagate) {
					//moving the old entry in case we couldnt propogate it to the next gate
					resultMaximalMappingQueue.push(queueEntry);
				}
			}
			else {
				//then none of the 2 qubits in the gate are mapped
					//we need to extend this mapping by finding possible pairs of qubits
					//so we will first find eligible controls, and for each control an eligible target
				bool propagated = false;
				for (auto j = copyArchitectureNodeDegrees.begin();j != copyArchitectureNodeDegrees.end();j++) {
					if (j->second >= circuitNodeDegrees[toSatisfy->control] && mappedArchQubits.find(j->second)==mappedArchQubits.end()) {
						//then this is a candidate for the control
						//we need to find adjacent qubits possible for the target
						auto innerRange = copyArchitectureEdges.equal_range(j->first);
						for (auto i = innerRange.first;i != innerRange.second;i++) {
							if (copyArchitectureNodeDegrees[i->second] >= circuitNodeDegrees[toSatisfy->target] && mappedArchQubits.find(i->second) == mappedArchQubits.end()) {
								//this pairing is a candidate for propagation
								propagated = true;
								map<int, int> newMapping = currMapping;
								set<int> newArchMapping = mappedArchQubits;
								set<GateNode*> newSatisfied = satisfiedGates;
								newMapping.insert(make_pair(toSatisfy->control, j->second));
								newMapping.insert(make_pair(toSatisfy->target, i->second));
								newArchMapping.insert(j->second);
								newArchMapping.insert(i->second);
								newSatisfied.insert(toSatisfy);
								resultMaximalMappingQueue.push(make_tuple(newMapping, newArchMapping, newSatisfied));
							}
						}
					}
				}
				if (!propagated) {
					//pushing old status then
					resultMaximalMappingQueue.push(queueEntry);
				}
			}

			
		}
		
		//setting our new queue to be the result queue (we basically just pushed propagations to a new queue and now we are setting it as the old queue)
		maximalMappingQueue = resultMaximalMappingQueue;
		otherClosedSet.insert(toSatisfy);
	}
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
	
	while (resultQueue.size() > 0) {
		tuple<map<int, int>, set<int>, set<GateNode*>> queueEntry = resultQueue.front();
		resultQueue.pop();
		if (get<2>(queueEntry).size() < mostMatched) {
			continue;
		}

		set<GateNode*> completedGates = get<2>(queueEntry);
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
bool wentTooFar(map<int,int>circuitNodeDegrees,map<int,int>architectureNodeDegrees,int circuitMaxDegree, int architectureMaxDegree) {
	//another prior check for going through both lists with two "fingers" and removing elements as I go for mappings
	//if I cannot remove all elements from the copyOfCircuitDegrees vector, then no complete initial mapping exists
	//we will need to roll back changes 
	//now we have left architecture nodes with at least the degree needed to satisfy the circuit requirements
	if (architectureNodeDegrees.size() < circuitNodeDegrees.size()) {
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
	for (auto i = architectureNodeDegrees.begin();i != architectureNodeDegrees.end();i++) {
		copyOfArchDegrees.push_back(i->second);
	}
	for (auto i = circuitNodeDegrees.begin();i != circuitNodeDegrees.end();i++) {
		copyOfCircuitDegrees.push_back(i->second);
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
		cout << "NO POSSIBLE MATCHING! SPOT 3\n";
		return true;
	}
	return false;
}

//adds an edge to our maximal mapping graph and returns the smallest and largest degree accordingly
pair<int,int> addGateEdge(GateNode* chosen, multimap<int,int>circuitEdges, map<int,int>circuitNodeDegrees) {
	//we add an edge from target to control, and from control to target in circuitEdges
	circuitEdges.insert(make_pair(chosen->control, chosen->target));
	circuitEdges.insert(make_pair(chosen->target, chosen->control));
	//incrementing degrees
	if (circuitNodeDegrees.find(chosen->control) != circuitNodeDegrees.end()) {
		circuitNodeDegrees[chosen->control]++;
	}
	else {
		circuitNodeDegrees[chosen->control] = 1;
	}
	if (circuitNodeDegrees.find(chosen->target) != circuitNodeDegrees.end()) {
		circuitNodeDegrees[chosen->target]++;
	}
	else {
		circuitNodeDegrees[chosen->target] = 1;
	}

	//returning the new smallest and greatest degrees of the resulting graph
	int maxDegree = 0;
	//setting min degree to greater than the max possible degree (number of vertices)
	int minDegree = circuitEdges.size();

	auto iter = circuitNodeDegrees.begin();
	while (iter != circuitNodeDegrees.end()) {
		if (iter->second < minDegree) {
			minDegree = iter->second;
		}
		if (iter->second > maxDegree) {
			maxDegree = iter->second;
		}

	}
	return std::make_pair(minDegree, maxDegree);
}

//given a max and min degree of the circuit graph I have defined, we can prune the coupling graph
	//returns all the changes made to the graph (in case we need to roll back)
 queue<pair<int,int>> pruneArchitectureEdges(multimap<int, int>architectureEdges, map<int, int>architectureNodeDegrees, int circuitMinDegree) {
	//iterating through architectureNodeDegrees and removing any vertex from the graph that is below the min degree	
	 //recording list of changes made
	 queue<pair<int, int>> changes;
	auto iter = architectureNodeDegrees.begin();
	while (iter != architectureNodeDegrees.end()) {
		if (iter->second < circuitMinDegree) {

			//then this node will not contribute to the maximal mapping
			
			//we will remove each edge one-by-one and as we remove, we will push our changes to the changes queue
			auto outerVertexIter = architectureEdges.equal_range(iter->first);
		    //we also need to remove the corresponding edge in the dst node
			auto edgeIter = outerVertexIter.first;
			while (edgeIter != outerVertexIter.second) {
				int dst = edgeIter->second;
				//firstly recording change to make
				//we only need to record one edge for each pair
				// this is because if we need to reverse changes, we can just remap (src,dest) and (dest,src) so we can save space
				changes.push(make_pair(iter->first, dst));
				auto innerRange = architectureEdges.equal_range(dst);
				auto innerIter = innerRange.first;
				while (innerIter != innerRange.second) {
					if (innerIter->second == edgeIter->first) {
						//reducing degree of this vertex
						architectureNodeDegrees[innerIter->second]--;
						if (architectureNodeDegrees[innerIter->second] == 0) {
							//then we can just remove the entire vertex also
							architectureEdges.erase(innerIter->second);
							architectureNodeDegrees.erase(innerIter->second);
							break;
						}
						//we found the edge to delete
						architectureEdges.erase(innerIter);
						
						break;
					}
					innerIter++;
				}
				
				//edgeIter=architectureEdges.erase(edgeIter);
			}
			//erasing entire vertex from degrees now and from graph

			architectureNodeDegrees.erase(iter->first);
			outerVertexIter = architectureEdges.equal_range(iter->first);
			architectureEdges.erase(outerVertexIter.first, outerVertexIter.second);
			//resetting everything since our degrees have been changed now
			iter = architectureNodeDegrees.begin();
		}
		else {
			
			iter++;
		}
	}
	//returning all the edges we erased
	return changes;
}

 //rolling back changes if we went too far in pruning
 void rollBack(queue<pair<int,int>>changes, multimap<int,int>architectureEdges, map<int,int>architectureNodeDegrees) {
	 while (changes.size() > 0) {
		 pair<int, int> item = changes.front();
		 changes.pop();
		 int src = item.first;
		 int dst = item.second;
		 //re-adding edges
		 architectureEdges.insert(make_pair(src, dst));
		 architectureEdges.insert(make_pair(dst, src));
		 //restoring degree
		 if (architectureNodeDegrees.find(src) != architectureNodeDegrees.end()) {
			 architectureNodeDegrees[src]++;
		 }
		 else {
			 architectureNodeDegrees[src] = 1;
		 }

		 if (architectureNodeDegrees.find(dst) != architectureNodeDegrees.end()) {
			 architectureNodeDegrees[dst]++;
		 }
		 else {
			 architectureNodeDegrees[dst] = 1;
		 }
	}

	//now all changes have been reverted
}

 //creating an initial array of distances of the shortest distance between any 2 physical qubits
 vector<vector<int>> getPhysicalDistancesFloydWarshall(multimap<int,int> architectureEdges) {
		//we run the vector based floyd warshall algorithm here for shortest distances	
	 //for compilation
	 vector<vector<int>> e;
	 return e;
}

 //naive distance implementation for now for distance between 2 physical qubits
 int getNaiveDistance(multimap<int, int> architectureEdges, int Q1, int Q2) {
		//we basically just use bfs and get shortest distance
	 queue<pair<int,int>> distanceQueue;
	 //initially putting first edges in queue
	 map<int,int> distances;
	 set<int> visited;
	 distances[Q1] = 0;
	 //int distance = 1;
	 auto range = architectureEdges.equal_range(Q1);
	 auto iter = range.first;
	 //int lastSource = Q1;
	 visited.insert(Q1);
	 while (iter != range.second) {
		 distanceQueue.push(make_pair(iter->first, iter->second));
		 iter++;
	 }
	 while (distanceQueue.size() > 0) {
		 pair<int, int> edge = distanceQueue.front();
		 distanceQueue.pop();
		 if (visited.find(edge.second) == visited.end()) {
			 //then this is the first time we are visiting this
			 //this is shortest distance
			 distances[edge.second] = distances[edge.first] + 1;
			 if (edge.second == Q2) {
				 return distances[edge.second];
			 }

			 //adding this to visited now
			 visited.insert(edge.second);

			 //adding neighbors to queue
			 range = architectureEdges.equal_range(edge.second);
			 iter = range.first;
			 while (iter != range.second) {
				 if (visited.find(iter->second) == visited.end()) {
					 //then we add to queue
					 distanceQueue.push(make_pair(iter->first, iter->second));
				 }
				 iter++;
			 }
		 }
	 }
	 //for compilation
	 return 0;
 }

 //used alongside floyd warshall algorithm to generate the distance matrix
 int getSmartDistance(vector<vector<int>> distanceMatrix, int Q1, int Q2) {
	 return distanceMatrix[Q1][Q2];
 }

 int getDistance(vector<vector<int>> distanceMatrix, multimap<int,int> architectureEdges, int Q1, int Q2) {
	 if (SMART) {
		 return getSmartDistance(distanceMatrix, Q1, Q2);
	 }
	 else {
		 return getNaiveDistance(architectureEdges, Q1, Q2);
	 }
 }



//helper to get my estimated partial mapping cost
 int partialMappingCost(map<int,int> mapping, set<GateNode*> remainingGates, vector<vector<int>> distanceMatrix, multimap<int,int> architectureEdges) {
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
		 }
		 else if (mapping.find(toAppraise->target) != mapping.end()) {
			 //then only target is mapped to some physical qubit
			 //we check if all neighbors of the target in the hardware mapping are mapped already
				//then we need to consider at least 1 swap in the cost (because there is no room to just assign the control mapping)
		 }

		 if (toAppraise->controlChild != NULL ) {
				//then we need to check if target parent is in the closedSet or not for insertion (so we enforce level order)
				if (closedSet.find((toAppraise->controlChild)->targetParent) != closedSet.end()) {
					//then we can add this child to the queue
					remainingDependencies.push(toAppraise->controlChild);
				}
			}

			if (toAppraise->targetChild != NULL ) {
				//then we need to check if control parent is in the closedSet or not for insertion (so we enforce level order)
				if (closedSet.find((toAppraise->targetChild)->controlParent) != closedSet.end()) {
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
 int partialMappingStitchingCost(map<int,int> firstMap, map<int,int> secondMap,vector<vector<int>> distanceMatrix, multimap<int,int> architectureEdges) {
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
	 }
	 //now we have the minimum # (lower bound) on swaps needed in order to transform a mapping into another mapping
	 return maxDistance;

 }

 //function to get a list of swaps needed to transform one mapping to another
 void stitchMappings(pair<map<int,int>,set<int>> mapOne, pair<map<int,int>,set<int>> mapTwo, multimap<int,int> architectureEdges) {
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
	 priority_queue<tuple<int,
		 queue<pair<int,int>>,
		 map<int,int>,
		 map<int,int>>> swapQueue;
	 map<int, int> firstMap = mapOne.first;
	 set<int> firstArchMap = mapOne.second;
	 map<int, int> secondMap = mapTwo.first;
	 set<int> secondArchMap = mapTwo.second;
	 while (true) {
		 if (swapQueue.size() == 0) {
			 // then we have to insert all initial swaps into the queue first with their respective current cost + heuristic cost
			 for (auto i = firstMap.begin();i != firstMap.end();i++) {
				 if (secondMap.find(i->first) != secondMap.end()) {
					 //then we can consider a swap for 1 distance
						//we check all neighbors in architecture in distance 1
						//if a logical qubit is mapped to a neighbor, that is a potential swap
						//if the size of the end mapping is > size of start mapping, we can set an unmapped qubit as a potential start neighbor
							//and then do the swap
					 auto edgeRange = architectureEdges.equal_range(i->first);
					 for (auto j = edgeRange.first;j != edgeRange.second;j++) {
						 //getting architecture nodes that are distance 1 away
						 //firstly seeing if these nodes are mapped in the start mapping
						 if (firstMap.find(j->second) != firstMap.end()) {
							 //then this is a possible swap to take
						 }
						 else {
							 //then this target node of the edge is not mapped
							 //we may need to adjust our initial mapping if the target mapping has more qubits than the start mapping
							 if (secondMap.size() > firstMap.size()) {
								 //now we can assign a mapping to the logical qubit 
							 }
						 }
					 }
				}
			 }
			 continue;
		 }
		 //getting the highest value operation
	 }
 }

 //lazy function to insert swaps (just tacks them on at the end of a maximal partition)
 void lazySwapInsert() {

 }

 //smarter function to insert swaps (considers delay and latency) 
 void smartSwapInsert() {

 }
