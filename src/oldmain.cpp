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

//signatures
bool satisfiesCircuit(int circuitNodeToMatch,int architectureNodeToMatch, multimap<int,int>& circuitEdges, multimap<int,int>& architectureEdges,vector<int>& mapping);
void addNeighbors(int circuitNodeToMatch, int architectureNodeToMatch,vector<pair <int, int>>& architectureNodeDegrees, vector<pair<int,int>>& circuitNodeDegrees, multimap<int, int>& circuitEdges, multimap<int, int>& architectureEdges, vector<int>& mapping,vector<int>& archMapping, list<pair<int, int>>& toAdd);
 
int main(int argc, char** argv) {
	char * qasmFileName = NULL;
	char * couplingMapFileName = NULL;
	int latency1 = 1;
	int latency2 = 1;
	int latencySwp = 1;
	
	//Parse command-line arguments:
	for(int iter = 1; iter < argc; iter++) {
		if(!strcmp(argv[iter], "-latency")) {
			latency1 = atoi(argv[++iter]);
			latency2 = atoi(argv[++iter]);
			latencySwp = atoi(argv[++iter]);
		} else if(!qasmFileName) {
			qasmFileName = argv[iter];
		} else if(!couplingMapFileName) {
			couplingMapFileName = argv[iter];
		} else {
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
	//new idea: lets try and introduce swaps when we can 


	//debug
	cout<<"PRINTING COUPLING\n";
	for (auto i=couplings.begin();i!=couplings.end();i++){
		printf("%d connected to %d\n",(*i).first,(*i).second);
	}


	//creating a multimap for quick lookup of architecture edges
	multimap<int, int> architectureEdges;
	//idea create a collection of node-> degree pairs
	//every entry initialized to (0,0)
	vector<pair <int, int>> architectureNodeDegrees(numPhysicalQubits,std::make_pair(0,0));
	//filling out the degrees from the couplings set
	for (auto i =couplings.begin() ;i!=couplings.end();i++) {
		// incrementing associated source and destination in architectureNodeDegrees
		int source = (*i).first;
		int dest = (*i).second;
		
		
		architectureNodeDegrees[source].first = source;
		architectureNodeDegrees[source].second ++;

		architectureNodeDegrees[dest].first = dest;
		architectureNodeDegrees[dest].second ++ ;


	}

	//debug
	cout<<"PRINTING ARCHITECTURE DEGREES:\n";
	for (auto i=architectureNodeDegrees.begin();i!=architectureNodeDegrees.end();i++){
		printf("%d has degree %d\n",(*i).first,(*i).second);
	}

	//now architectureNodeDegrees is a vector of the form [(Q_0,DEGREE),(Q_1,DEGREE),(Q_2,DEGREE),...)]
	// also creating an unordered_set of (src,dest) node pairs for each gate we encounter in the circuit 
	multimap<int,int> circuitEdges;

	//So, for circuit creating node->degree pairs below (MAKE SURE TO KEEP TRACK OF MIN DEGREE)
	vector<pair<int, int>> circuitNodeDegrees(numLogicalQubits , std::make_pair(0, 0));
	auto iter = firstGates.begin();
	//creating unordered set for closed set checking
	unordered_set<GateNode*> closedSet;
	while (iter != firstGates.end()) {
		//creating stack for dfs
		stack<GateNode*> fringe;
		
		GateNode* toExplore = *iter;
		fringe.push(toExplore);

		while (!fringe.empty()) {
			GateNode* curr = fringe.top();
			fringe.pop();
			if (closedSet.count(curr)){
				//then we skip this 
				continue;
			}
			// getting control and target logical qubits
			int control = curr->control;
			int target = curr->target;
			if (control != -1) {
				//debug
				printf("control: %d and target: %d\n",control,target);
				(circuitNodeDegrees[control]).first = control;
				//increasing degree
				(circuitNodeDegrees[control]).second ++;

				//target is always present 
				(circuitNodeDegrees[target]).first = target;
				(circuitNodeDegrees[target]).second ++;

				//adding this control target pair to our unordered_map of edges in a bidirectional fashion
				circuitEdges.insert(make_pair(control, target));
				circuitEdges.insert(make_pair(target, control));

				
			}
			// if control==-1, then this is just a single qubit gate, which can execute on itself	
			// incrementing their degree values in our vector
			GateNode* targetChildren = curr->targetChild;
			GateNode* controlChildren = curr->controlChild;

			if (targetChildren != NULL) {
				if (closedSet.count(targetChildren)==0) {
					// then I add to the fringe
					fringe.push(targetChildren);
				}
			}

			if (controlChildren != NULL) {
				if (closedSet.count(controlChildren) == 0) {
					// then I add to the fringe
					fringe.push(controlChildren);
				}
			}
						

			//adding the explored node to the closed set
			closedSet.insert(curr);
		
		}

		//incrementing iterator
		iter++;
	}

	//now we have logical qubits and their degrees 
	
	//debug
	cout<<"Printing circuit degrees:\n";
	for (auto i=circuitNodeDegrees.begin();i!=circuitNodeDegrees.end();i++){
		printf("%d has degree %d\n",(*i).first,(*i).second);
	}

	//debug
	cout<<"PRINTING CIRCUIT EDGES\n";
	for (auto i=circuitEdges.begin();i!=circuitEdges.end();i++){
		printf("%d connected to %d\n",(*i).first,(*i).second);
	}

	int minDegreeEncountered;
	// getting smallest degree in the vector
	for (int i = 0;i < numLogicalQubits;i++) {
		if (i == 0) {
			minDegreeEncountered = circuitNodeDegrees[i].second;
			continue;
		}
		int currDegree = circuitNodeDegrees[i].second;
		if (currDegree < minDegreeEncountered) {
			minDegreeEncountered = currDegree;
		}
	}
	


	//Removing nodes from architecture with a degree less than min degree of any node in the circuit, this will help in computation later
	for (auto i = architectureNodeDegrees.begin();i != architectureNodeDegrees.end();i++) {
		pair<int, int> toConsider = *i;
		if (toConsider.second < minDegreeEncountered) {
			// then we remove this
			
			i = architectureNodeDegrees.erase(i);

			//need to adjust degree of neighbors to be 1 less now after removal


		}
	}

	//intializing map of architecture edges now (after pruning nodes that will not contribute to our circuit)
	for (auto i = architectureNodeDegrees.begin();i != architectureNodeDegrees.end();i++) {
		int nodeToAdd = (*i).first;
		for (auto j = couplings.begin();j != couplings.end();j++) {
			pair<int, int> srcDest = (*j);
			int src = (*j).first;
			int dst = (*j).second;
			if (src == nodeToAdd ) {
				architectureEdges.insert(make_pair(src, dst));
			}
			else if (dst== nodeToAdd) {
				architectureEdges.insert(make_pair(dst, src));
			}
		}
	}

	//debug
	cout<<"PRINTING ARCHITECTURE EDGES\n";
	for (auto i=architectureEdges.begin();i!=architectureEdges.end();i++){
		printf("%d connected to %d\n",(*i).first,(*i).second);
	}

	//now we have left architecture nodes with at least the degree needed to satisfy the circuit requirements
	if (architectureNodeDegrees.size() < circuitNodeDegrees.size()) {
		//then we do not have the required number of nodes to try and find a perfect initial mapping
		cout << "NO MAPPING POSSIBLE! SPOT 1\n"  ;
		return 0;
	}

	//sorting copies of both lists from least to greatest node degrees
	vector<pair<int, int>>copyOfArchDegrees = architectureNodeDegrees;
	vector<pair<int, int>>copyOfCircuitDegrees = circuitNodeDegrees;
	sort(copyOfArchDegrees.begin(), copyOfArchDegrees.end());
	sort(copyOfCircuitDegrees.begin(), copyOfCircuitDegrees.end(),
		[](pair<int, int>a, pair<int, int>b) {return a.second > b.second;});


	//another prior check for going through both lists with two "fingers" and removing elements as I go for mappings
	//if I cannot remove all elements from the copyOfCircuitDegrees vector, then no complete initial mapping exists
	auto architectureNodesFinger = copyOfArchDegrees.begin();
	auto circuitNodesFinger = copyOfCircuitDegrees.begin();

	while (architectureNodesFinger != copyOfArchDegrees.end() &&
		circuitNodesFinger != copyOfCircuitDegrees.end()) {
		if ((*circuitNodesFinger).second <= (*architectureNodesFinger).second) {
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
		cout << "NO POSSIBLE MATCHING! SPOT 2\n";
		return 0;
	}

	//see if this actually does something or not
	//copyOfArchDegrees.clear();
	//copyOfCircuitDegrees.clear();


	//need to maintain a listing of neighbors required for each node in both the architecture and in the circuit
	
	//initializing a mapping vector
	//-1 means this is not mapped yet
	vector<int> mapping(numLogicalQubits, -1);
	vector<int> archMapping(numPhysicalQubits, -1);
	

	bool done = false;


	vector<pair<int, int>> lastMapping;

	//keeping track of the best mapping during this process (most nodes mapped with no issues)
	//this will help in getting an initial mapping for the actual project
	int mostMapped = 0;
	vector<int> bestMappingSoFar;

	
	
	//other idea: at each step, I have a list of pairs to match if ever a list is empty in a spot, we backtrack to the last spot
	//so we have something like this:
	//this is basically a circularly linked list, where the data is a reference to another circularly linked list
		//this is for matching
	list<pair<int, int>> matchingStructure;



	//keeping track of our levels
	vector<std::list<pair<int,int>>::iterator> levelIterators;
	
	int smallestDegreeNode = copyOfCircuitDegrees[0].first;

	auto archIterator = architectureNodeDegrees.begin();
	while (archIterator != architectureNodeDegrees.end()) {
	
		matchingStructure.push_front(make_pair(smallestDegreeNode, (*archIterator).first));
		archIterator++;
	}

	
	//very greedy idea below, kind of optimized it by swapping elements across levels instead of just repushing it

	while (!matchingStructure.empty()) {
		//debug
		cout<<"PRINTING STACK\n";
		for (auto i=matchingStructure.begin();i!=matchingStructure.end();i++){
			printf("stack: (%d,%d)\n",(*i).first,(*i).second);
		}


		//getting elements		

		
		pair<int, int> suggestedMatching = matchingStructure.front();
		//removing the pair from our actual stack
		matchingStructure.pop_front();
		
		if ( levelIterators.size()!=0 && (*(levelIterators[levelIterators.size() - 1])) == suggestedMatching) {
			printf("WE ARE BACKTRACKING!\n");
			if (lastMapping.size()> mostMapped) {
				bestMappingSoFar = mapping;
				mostMapped=lastMapping.size();
			}
			else if (lastMapping.size() == mostMapped) {
				//determining number of swaps needed to satisfy the circuit and picking best mapping so far
			}
			//then we are backtracking and nothing in the next level worked
			//so we rollback a matching
			// keep track of best matching so far
			pair<int, int> lastMatched = lastMapping.back();
			int src = lastMatched.first;
			int dst = lastMatched.second;
			mapping[src] = -1;
			archMapping[dst] = -1;
			//erasing the last element 
			lastMapping.pop_back();
			//erasing the iterator
			levelIterators.pop_back();
			continue;
		}
		
		int circNodeToMatch = suggestedMatching.first;
		int archNodeToMatch = suggestedMatching.second;

	

		if (circuitNodeDegrees[circNodeToMatch].second > architectureNodeDegrees[archNodeToMatch].second) {
			//then we cannot match these at all
			continue;
		} 

		//checking if matching is satisfied if we reach here

		//mapping
		mapping[circNodeToMatch] = archNodeToMatch;
		archMapping[archNodeToMatch] = circNodeToMatch;

		if (!satisfiesCircuit(circNodeToMatch,archNodeToMatch,circuitEdges,architectureEdges,mapping)) {
			//then we can remove the mapping and continue because circuit specifications were violated
			mapping[circNodeToMatch] = -1;
			archMapping[archNodeToMatch] = -1;
			continue;
				
		}
		printf("Mapped %d to %d \n",circNodeToMatch,archNodeToMatch);
		//if we reached here, then mapping is good so far

		//putting successful mapping onto list of mappings
		lastMapping.push_back(make_pair(circNodeToMatch, archNodeToMatch));

		//getting iterator to this element to push
		std::list<pair<int, int>>::iterator toAdd = matchingStructure.begin();

		//greedy part of solution that needs to be optimized
		//idea for optimization: repush the nodes, but also remove them from the bottom because they represent a different pathway
		if (levelIterators.size()!= 0) {
			auto iterator = levelIterators.back();
			printf("back of iterators: (%d,%d)\n",(*iterator).first,(*iterator).second);
			while (iterator != toAdd) {
				pair<int, int> toConsider = (*iterator);
				//if the src or dst has already been matched, then we do not repush it
				if (mapping[toConsider.first] == -1 && archMapping[toConsider.second] == -1) {
					//then we push this
					matchingStructure.push_front(toConsider);
					printf("REPUSHING (%d,%d)\n",(toConsider).first,(toConsider).second);
					//deleting logic
					
					iterator=matchingStructure.erase(iterator);
					continue;
					
				}
				iterator--;
			}
				
		}
		
		

		//before adding neighbors we denote this iterator as the last level
		levelIterators.push_back(toAdd);
		printf("Iterator placed at (%d,%d)\n",(*toAdd).first,(*toAdd).second);
	
		

		//adding neighbors of circNode that was matched and archNode that was matched
		addNeighbors(circNodeToMatch, archNodeToMatch,circuitNodeDegrees,architectureNodeDegrees, circuitEdges, architectureEdges, mapping, archMapping, matchingStructure);

		cout<<"Printing current matching formed!\n";
		for (int i=0;i<lastMapping.size();i++){
			printf("%d matched to %d\n",lastMapping[i].first,lastMapping[i].second);
		}

		if (lastMapping.size() == numLogicalQubits) {
			//then we are done with the mapping
			done = true;
			break;
		}


	}

	if (done) {
		//then we have successful mapping that we can print out

		printf("Found a perfect initial mapping!\n");
		for (int i = 0;i < mapping.size();i++) {
			printf("q%d mapped to Q%d\n", i, mapping[i]);
		}
	}
	else {
		cout << "NO PERFECT INITIAL MAPPING POSSIBLE!\n";
		if (bestMappingSoFar.size()){
			cout<< "Here is the best mapping we found:\n";
			for (int i=0;i<bestMappingSoFar.size();i++){
				if (bestMappingSoFar[i]!=-1){
					//then this is mapped
					printf("q%d mapped to Q%d\n", i, bestMappingSoFar[i]);
				} else{
					printf("q%d is unmapped!\n",i);
				}
			}
		}
	}
	

	
	
	//Exit the program:
	return 0;
}


//helper function to add (circNode,archNode) pairs of neighbors of given nodes to our list
void addNeighbors(int circuitNodeToMatch, int architectureNodeToMatch,vector<pair <int, int>>& architectureNodeDegrees, vector<pair<int,int>>& circuitNodeDegrees, multimap<int, int>& circuitEdges, multimap<int, int>& architectureEdges, vector<int>& mapping,vector<int>& archMapping, list<pair<int, int>>& toAdd) {
	

		//pushing architectureNodeNeighbors

		
		auto circRange = circuitEdges.equal_range(circuitNodeToMatch);
		auto circIterator = circRange.first;
		while (circIterator != circRange.second) {
			if (mapping[(*circIterator).second] == -1) {
				// then we push
				auto archRange = architectureEdges.equal_range(architectureNodeToMatch);
				auto archIterator = archRange.first;
				while (archIterator != archRange.second) {
					if (archMapping[(*archIterator).second] == -1) {
						//then we confirm the push
						//CHECKING IF DEGREES ARE COMPATIBLE BEFORE ADDING THE PAIR
						if (circuitNodeDegrees[(*circIterator).second].second <= architectureNodeDegrees[(*archIterator).second].second) {	
							toAdd.push_front(make_pair((*circIterator).second,(*archIterator).second));
						}
					}
					archIterator++;
				}
			}
			circIterator++;
		}

		
}


//helper function for deciding if the matching is satisfied so far
bool satisfiesCircuit(int circuitNodeToMatch,int architectureNodeToMatch, multimap<int,int>& circuitEdges, multimap<int,int>& architectureEdges,vector<int>& mapping) {
	auto range = circuitEdges.equal_range(circuitNodeToMatch);
	auto circIter = range.first;
	while (circIter != range.second) {
		int otherMappedArchNode = mapping[(*circIter).second];
			if (otherMappedArchNode != -1) {
				// then we should check the neighbor of this node
				auto archRange = architectureEdges.equal_range(otherMappedArchNode);
				auto archIter = archRange.first;
				bool innerCheck = false;
					while (archIter != archRange.second) {
						//then we should check if mapping[i] is here
						if ((*archIter).second == architectureNodeToMatch) {
							//then this is satisfied and we move on to check next one
							innerCheck = true;
							break;
						}
						
						archIter++;

					}
					if (!innerCheck) {
						//then we failed
						return false;
					}
			}
		circIter++;
	}	
	// if we reached here, then the circuit is satisfied so far
	return true;

			

}
