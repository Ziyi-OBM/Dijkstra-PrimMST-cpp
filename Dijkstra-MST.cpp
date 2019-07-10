

/*

06/28/2019
This code append the code from homework 2 - Dijkstra with the following
1.Constructor to the Graph class that takes a istream and read the file into a graph of matrix representation. L:171~196
2.Method in the Graph class that finds the minimum spanning tree. L:593~669
3.New helper class called "spanTree" to wrap up the result of the MST algorithm. L:53~61, 86~106

In addition, a bug from hw2 submission was found and fixed:
*Fixed the problem that the shortest path algorithm only works when starting at node 0 (was hard coded this way in the algorithm)
*/


#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <fstream>
#include <iterator>
#include<numeric>


//A class used to wrap the information of a edge, i.e. starting node, ending node and cost.		
class edgeQ{
	public:
	
	edgeQ():cost(0.0),start(0),ends(0){}
	edgeQ(double cost, int start, int ends):cost(cost),start(start),ends(ends){}
	
	double cost=0.0;
	int start=0;
	int ends=0;
	
	~edgeQ(){}
};


//A class used to wrap the information of a shortest path solution, i.e. total cost and the predecessors
class pathData{
	public:
	
	double cost=0.0;
	std::vector<int> route;
	
	pathData():cost(-1.0),route(){}

	~pathData(){}
};
//A class used to wrap the information of a MST solution, i.e. total cost and the edges
class spanTree{
	public:
	
	double cost=0.0;	//The total cost of the MST
	std::vector<int> edge;	//Vector of size "# of nodes". Edges in the MST can be represented as {edge[n],n} for 0<=n<=number of vertices
	
	spanTree():cost(-1.0),edge(){}

	~spanTree(){}
};

std::ostream& operator<< (std::ostream& out, const edgeQ& edge){
	out << "Edge Object:\nCost: "<<edge.cost<<"; Start at: "<< edge.start<<"; Ends at: "<<edge.ends;
	return out;}
	
	
std::ostream& operator<< (std::ostream& out, const pathData& path){
	if (path.cost!=-1){
		std::cout << "The cost of the shortest path is: "<< path.cost <<std::endl;
		
		std::cout << "The route is: "  << std::endl;
		std::vector<int> toPrint = path.route;
		
		//https://stackoverflow.com/questions/18446946/printing-vector-in-reverse-order/54448537
		for(auto it = toPrint.crbegin(); it != toPrint.crend(); it++){
			std::cout<<" -> "<<*it;
		}
		std::cout << std::endl;
	}else if(path.cost==-1)
		std::cout << "No path exists."  << std::endl;
	return out;
}

std::ostream& operator<< (std::ostream& out, const spanTree& tree){
	if (tree.cost!=-1){
		std::cout << "The cost of the minimum spanning tree is: "<< tree.cost <<std::endl;
		
		std::cout << "The edges to reach each vertices are: "  << std::endl;
		std::vector<int> toPrint = tree.edge;
		
		//https://stackoverflow.com/questions/18446946/printing-vector-in-reverse-order/54448537
		int nodeNum=0;
		for(auto it:toPrint){
			if (it == nodeNum){
				std::cout<<"Vertex #"<< nodeNum <<" <--- The starting point of the tree" <<std::endl;
			}else{
				std::cout<<"Vertex #"<< nodeNum <<" <--- From Vertex #"<< it <<std::endl;
			}
			nodeNum++;
		}
	}else if(tree.cost==-1)
		std::cout << "Graph is not connected. Not MST"  << std::endl;
	return out;
}


class Graph{
		

	inline double prob(){return (rand() / ((double)RAND_MAX+1));}
	template<typename T>
	inline void disp(T vec2D) const {	//Print vector borrowed from stack overflow
		for(auto vec : vec2D){
			for(auto x : vec)
				std::cout<<x<<" , ";
			std::cout << std::endl;
		}
	}
	template<typename T>
	inline int count(T vec2D) const {
		int sum = 0;
		for(auto vec : vec2D)
		{
			for(auto x : vec)
				if(x!=0.0)++sum;
		}
		return sum;
	}

	int V;	//# of vertices	
	//Initialize a VxV bool matrix to represent the nodes
	std::vector<bool> rowg;
	std::vector<std::vector<bool> >   G;

	//Initialize a VxV double matrix to store the cost of the edges
	std::vector<double> rowe;
	std::vector<std::vector<double> >  E;
	
	public:
	//Parameters:
	
	
	//Default constructor
	Graph():V(1),rowg(V),G(V,rowg),rowe(V,0.0),E(V,rowe){
		
		std::cout << "Default constructor called" << std::endl;
		
		// //Initialize a VxV bool matrix to represent the nodes
		// rowg = std::vector<bool>(V);
		// G = std::vector< std::vector<bool> >(V,rowg);

		// //Initialize a VxV double matrix to store the cost of the edges
		// rowe = std::vector<double>(V,0.0);
		// E = std::vector< std::vector<double> >(V,rowe);		
	}

	//Parameterized constructor
	Graph(int ver):V(ver),rowg(V),G(V,rowg),rowe(V,0.0),E(V,rowe){
		
		// //Initialize a VxV bool matrix to represent the nodes
		// rowg = std::vector<bool>(V);
		// G = std::vector< std::vector<bool> >(V,rowg);

		// //Initialize a VxV double matrix to store the cost of the edges
		// rowe = std::vector<double>(V,0.0);
		// E = std::vector< std::vector<double> >(V,rowe);		
	}
	
	//Read File cnstructor	
	Graph(std::ifstream *inFile){	
		std::cout << "Read constructor called"  << std::endl;
		std::vector<int> inGraph;
		int val;
		while (*inFile >> val)
		{
			inGraph.push_back(val);
		}
		
		V=inGraph[0];
		int edgeNum=inGraph.size()-1;
		
		//Initialize a VxV bool matrix to represent the nodes
		rowg = std::vector<bool>(V);
		G = std::vector< std::vector<bool> >(V,rowg);

		//Initialize a VxV double matrix to store the cost of the edges
		rowe = std::vector<double>(V,0.0);
		E = std::vector< std::vector<double> >(V,rowe);	

		for (int i=1 ; i < edgeNum-3; i=i+3){
			addEdge (inGraph[i], inGraph[i+1], inGraph[i+2]);			
		}		
		
	}

	
	

	friend std::ostream& operator<< (std::ostream& out, const Graph& g);
	
	//Function to randomly initialize the graph 
	void generateRanomGraph(double edgeDensity, double distRangeLb, double distRangeUb) { 
		std::cout << "Randomly generating graph of density: " << edgeDensity << " and distance range: " << distRangeLb <<" ~ "<< distRangeUb <<std::endl; 
		srand(clock());   
		// srand(20190612);  Setting a constant seed for debugging
		if (distRangeLb>distRangeUb){
			std::cout << "Invalidd Distance Bounds" << std::endl;
			return;
		}
		double range = distRangeUb - distRangeLb;
		  for (int i=0;i<V;++i){
			for (int j=i;j<V;++j){
					if(i==j)G[i][j]=false;
						else{
							G[i][j]=G[j][i]=(prob()<edgeDensity);
							//if(G[i][j])E[i][j]=E[j][i]=(distRange*(1-prob()));
							if(G[i][j])E[i][j]=E[j][i]=(range*prob()+distRangeLb);
						}
			}
		  }	  
		
		//Print the resulting array to debug
		//disp(G);
		//std::cout << "\n" << std::endl;
		//disp(E);	   
	} 

	//Interface methods
	void printNodes() const {disp(G);}	//Print the Nodes as a matrix

	void printEdges() const {disp(E);}	//Print the Edges as a matrix

	int vertices() const {return V;} //Returns the number of vertices in the graph

	int edges() const {return count(E);}	//Returns the number of edges in the graph

	bool adjacent(int i, int j) const {return G[i][j];}	//Tests whether there is an edge from node i to node j.

	std::vector<int> neighbors(int i) const {	//Print all nodes j such that there is an edge from i to j.
		std::vector<int> neighborNodes;
		for (int j=0;j<V;++j){
			if(G[i][j]){
				neighborNodes.push_back(j);
			}
		}
		return neighborNodes;
	}

	void addEdge (int i, int j, double cost){	//adds to G the edge from i to j with cost x, if it is not there.
			if(!G[i][j]){
				G[i][j]=G[j][i]=true;
				E[i][j]=E[j][i]=cost;
			}//else{std::cout<<"Error: Edge exists. The value is unchanged. Use setEdgeValue to assign new value"<<std::endl;}
	}

	void deleteEdge (int i, int j){	//removes the edge from i to j, if it exists.
			if(G[i][j]){
				G[i][j]=G[j][i]=false;
				E[i][j]=E[j][i]=0.0;
			}else{std::cout<<"Warning: Edge does not exists. No deletion peformed"<<std::endl;}
	}

	bool getNodeValue(int i, int j) const {return G[i][j];}	//returns the value associated to the node (i,j) (True/False=there Is/Not a node).

	// set_node_value is not necessary because its the same as the "addEdge" method

	double getEdgeValue(int i, int j) const {return E[i][j];}	//returns the value associated to the edge (i,j).

	double setEdgeValue(int i, int j, double cost){	//sets the value associated to the edge (i,j) to v.
		//If the two nodes are not connected, remind the user to use "addEdge" to add a new connection between the two nodes
		if(!G[i][j])std::cout<<"Error: No edge exists between given nodes. please use the 'addEdge' method"<<std::endl;
						else{E[i][j]=E[j][i]=cost;}
	}

	spanTree primMST(int u=0);
			
		
	//Destructor
	~Graph() {} 


	//Private variables
	private:




		
};

std::ostream& operator<< (std::ostream& out, const Graph& g){
	std::cout << "Graph Object: "  << std::endl;
	std::cout << "Node Matrix: "  << std::endl;
	for(auto vec : g.G){
			for(auto x : vec)
				std::cout<<x<<" , ";
			std::cout << std::endl;
		}
		
	std::cout << "Edge Matrix: "  << std::endl;
	for(auto vec : g.E){
			for(auto x : vec)
				std::cout<<x<<" , ";
			std::cout << std::endl;
		}	
		
	return out;
}
	
	
	
	
class priorityQueue{
	//Implement priority queue using min heap	
	//Reference: https://www.geeksforgeeks.org/binary-heap/
	template<typename T>
	inline void swap(T *a, T *b){	//Print vector borrowed from stack overflow
		T temp = *a;
		*a=*b;
		*b=temp;
	}	

	inline int parentId(int i) { return (i-1)/2; } 
	inline int LchildId(int i) { return (2*i)+1; } 
	inline int RchildId(int i) { return (2*i)+2; } 
	inline int smaller(int i,int j) const {
		if (minHeap[i].cost<minHeap[j].cost){
			return i;
		}else{
			return j;
		}
		
	}

	//Defining member variables
	int capacity;
	edgeQ edgeInit=edgeQ();
	std::vector<edgeQ> minHeap;
	int heapEnd=0;	
	
	public: 	
	//Default constructor: Create an empty queue of max size 5000
	priorityQueue():capacity(5000), minHeap(capacity,edgeInit) {}
	//Constructor: Create an empty queue of max size cap
	priorityQueue(int cap):capacity(cap), minHeap(capacity,edgeInit) {}
	
	void getCapacity() const { // Print capacity
				std::cout << "The capacity of the Queue is: " << capacity <<std::endl; 
	}
	
	//Test if the queue contains a path to a certain destination
	//Returns the index of the edge in the queue, or return -1 if it does not contain such path
	int contains(int destination) const{
		int itIndex=0;
		for (auto it : minHeap){
			if (it.ends == destination){
				return itIndex;
			}
			itIndex++;
		}
		return -1;
	}
	
	//Increase the priority(decrease the cost) of a edge in the queue which has index: nodeIndex in the minHeap vector, 
	//also update the start point of this new edge/path 
	void increasePriority (int nodeIndex, int newStart, double newCost){
		if (nodeIndex==-1){
			return;
		}
		if (minHeap[nodeIndex].cost<newCost){
			//std::cout<<"The updated route does not improve from the existing one."<<std::endl;
			return;
		}		
		minHeap[nodeIndex].start=newStart;
		minHeap[nodeIndex].cost=newCost;		
		moveUp(nodeIndex);
		return;
	}
		
	
	void insert(edgeQ edgeIn){ // insert queue_element into queue
		if (heapEnd == capacity){ 
				std::cout << "\nInsertion failed: Queue reached max capacity\n" << std::endl; 
				return;
			} 
		minHeap[heapEnd]=edgeIn;	

		moveUp(heapEnd);
		heapEnd++; 
	}
	
	edgeQ extract(){ // Return and delete the root element
		if (heapEnd == 0){ 
				std::cout << "\nExtract failed: Queue is empty\n" << std::endl; 
				return edgeQ();
			} 
			
		if (heapEnd == 1){ 
				heapEnd--;
				return minHeap[0];
			} 

		edgeQ root = minHeap[0];

		
		heapEnd--;
		minHeap[0]=minHeap[heapEnd];

		
		moveDown(0);
		
		return root;
	}

	edgeQ top() const { // Return the root element
		if (heapEnd == 0){ 
				std::cout << "\nExtract failed: Queue is empty\n" << std::endl; 
				return edgeQ();
			} 
		
		return minHeap[0];
	}
	
	
	int sizeQ() const { // Return the root element
		std::cout << "The size is: " << heapEnd << "\n" << std::endl; 
		return heapEnd;
	}
	
	bool isEmpty() const { // Return the root element
		if (heapEnd == 0){ 
			return true;
		}else{
			return false;
		}
	}
	
	void clearAll(){
		heapEnd=0;		
	}

	//Destructor
	~priorityQueue() {} 


	private:

	//Compare a node with its parents and move it upward to maintain heap property
	void moveUp(int cursor){
	int parent = parentId(cursor);
		while(cursor!=0 & minHeap[parent].cost>minHeap[cursor].cost){
			
			swap(&minHeap[parent],&minHeap[cursor]);
			
			cursor = parent;
			parent = parentId(cursor);

		}
	}
	
	//Compare a node with its childs and move it downward to maintain heap property
	void moveDown(int cursor){
		
		if (LchildId(cursor)>=heapEnd){
			return;
		}
		
		int smallerChild= LchildId(cursor);
		
		if (RchildId(cursor) < heapEnd & minHeap[smallerChild].cost > minHeap[RchildId(cursor)].cost){
			smallerChild = RchildId(cursor);
		}
		
		if ( minHeap[cursor].cost > minHeap[smallerChild].cost){
			swap(&minHeap[cursor],&minHeap[smallerChild]);	
			
			cursor=smallerChild;
			moveDown(cursor);
		}
		
	}

};


//Class to find the shortest path between two verties, given a graph
class shortestPath{
	
	Graph const *G;
	
	int vertices;
	public:
	
	//Initialize the object
	shortestPath(Graph const *graphIn):G(graphIn),vertices(G->vertices()){
		//Debugging: Printing the input information
		//std::cout << "Number of Nodes in the Graph: " << vertices <<std::endl;
		//std::cout << "Number of Edges in the Graph: " << G->edges() <<std::endl;
		//Q.getCapacity();
		//std::cout << *G <<std::endl;
	}
	
	//Method to find the shortest path
	pathData findPath(int u, int w){
		//Create a priority queue that serves as the openset discussed in the algorithm
		priorityQueue Q(G->edges());
		
		
		std::vector<int> neighborNodes = G->neighbors(u);
		std::vector<double> dist(vertices,-1);	//A vector that holds the shortest distance from the start point to each node
		std::vector<int> prev(vertices,-1);	//A vector that holds the previous nodes on the shortest path
		
		int closedSetSize = 0; 	//The close set is the visited nodes that we've found the shorted distance
								//It is the number of non-"-1" element in the dist and prev vector
		dist[u]=0;
		prev[u]=u;	
		closedSetSize++;
		
		//Add the neibours of the starting node to the openset
		for (auto x : neighborNodes){
			edgeQ edgeIn=edgeQ(G->getEdgeValue(u,x),u,x);
			if (dist[edgeIn.ends]==-1)
				Q.insert(edgeIn);
		//Debugging: Printing elements that get added to the openset.
		//std::cout << "Inserting: " << std::endl;
		//std::cout << edgeIn << std::endl;			
			
		}
		
		//While we have not reach destination and the openset is not exhausted, carry on the algorithm
		while (dist[w] == -1 & !Q.isEmpty()){
			
			//Extract the highest priority element in the openset and move it to the closed set
			edgeQ currentEdge=Q.extract();
			//Debugging: Printing elements that get added to the close set.
			//std::cout << "Evaluating: " << std::endl;
			//std::cout << currentEdge << std::endl;
		
			//Save the distance and path to the dist and prev vector
			dist[currentEdge.ends]=currentEdge.cost;
			closedSetSize++;
			prev[currentEdge.ends]=currentEdge.start;
			
			//Explore the neibours of the point that was added to the closed set
			neighborNodes = G->neighbors(currentEdge.ends);
			for (auto x : neighborNodes){
				//Gather the information of the edges to its neibour
				edgeQ edgeIn(G->getEdgeValue(currentEdge.ends,x)+dist[currentEdge.ends],currentEdge.ends,x);
				
				//Update the openset base on these neibours
				int isInQ = Q.contains(edgeIn.ends);
				if (isInQ==-1 & dist[edgeIn.ends]==-1) Q.insert(edgeIn);			
				//This method will test whether it improves over the previous path, and update it when it does			
				else if (isInQ!=-1 & dist[edgeIn.ends]==-1) Q.increasePriority (isInQ, edgeIn.start, edgeIn.cost);
				
			}			
		}
		
		std::vector<int> route;	//A variable to hold the vertices along the shortest path
		
		//After the growing process
		//First test if we have reached the destination
		if (dist[w]!=-1){
			//If true, wrap the result to a helper class
			int cursor = w;
			route.push_back(w);
			while (cursor != u){
				route.push_back(prev[cursor]);
				cursor=prev[cursor];
			}
			double totalCost = dist[w];	
			pathData pathOut;
			pathOut.route = route;
			pathOut.cost = totalCost;
			return pathOut;
			//If destination is not reached, return a empty path with cost -1
		}else if (Q.isEmpty()) return pathData();
		
	}	
	
	//Destructor
	~shortestPath() {} 
		
	private:

	
	
};


	//Implement a method to calculate the MST using prim's algorithm
	spanTree Graph::primMST(int u){
		
		int vertices = V;
		//Create a priority queue that serves as the openset discussed in the algorithm
		priorityQueue Q(edges());
		
		
		std::vector<int> neighborNodes = neighbors(u);
		std::vector<double> dist(vertices,-1);	//A vector that holds the shortest distance from the start point to each node
		std::vector<int> prev(vertices,-1);	//A vector that holds the previous nodes on the shortest path
		
		int closedSetSize = 0; 	//The close set is the visited nodes that we've found the shorted distance
								//It is the number of non-"-1" element in the dist and prev vector
		dist[u]=0;
		prev[u]=u;	
		closedSetSize++;
		
		//Add the neibours of the starting node to the openset
		for (auto x : neighborNodes){
			edgeQ edgeIn=edgeQ(getEdgeValue(u,x),u,x);
			if (dist[edgeIn.ends]==-1)
				Q.insert(edgeIn);
		//Debugging: Printing elements that get added to the openset.
		//std::cout << "Inserting: " << std::endl;
		//std::cout << edgeIn << std::endl;			
			
		}
		
		//While we have not reach all the nodes and the openset is not exhausted, carry on the algorithm
		while (closedSetSize<vertices & !Q.isEmpty()){
			
			//Extract the highest priority element in the openset and move it to the closed set
			edgeQ currentEdge=Q.extract();
			//Debugging: Printing elements that get added to the close set.
			//std::cout << "Evaluating: " << std::endl;
			//std::cout << currentEdge << std::endl;
		
			//Save the distance and path to the dist and prev vector
			dist[currentEdge.ends]=currentEdge.cost;
			closedSetSize++;
			prev[currentEdge.ends]=currentEdge.start;
			
			//Explore the neibours of the point that was added to the closed set
			neighborNodes = neighbors(currentEdge.ends);
			for (auto x : neighborNodes){
				//Gather the information of the edges to its neibour
				
				//Here is the essential difference between prim and Dijkstra, where in Dijkstra, 
				//the first argument, which is the cost to the node would be getEdgeValue(currentEdge.ends,x)+dist[currentEdge.ends]
				//Which is the cost from the source to x. Here instead, the cost is the from the point we are exploring to x.
				edgeQ edgeIn(getEdgeValue(currentEdge.ends,x),currentEdge.ends,x);
				
				//Update the openset base on these neibours
				int isInQ = Q.contains(edgeIn.ends);	//Return the node number if contains and return -1 if not contain
				//If not in closed set and NOT in open set, insert the edge
				if (isInQ==-1 & dist[edgeIn.ends]==-1) Q.insert(edgeIn);		
				//If not in closed set but IS in open set, compare it and replace the edge if this edge has a lower cost	
				else if (isInQ!=-1 & dist[edgeIn.ends]==-1) Q.increasePriority (isInQ, edgeIn.start, edgeIn.cost);
				
			}			
		}
		
		std::vector<int> tree;	//A variable to hold the vertices along the shortest path
		
		//After the growing process
		//First test if we have reached all the vertices
		if (closedSetSize==vertices){
			//If true, wrap the result to a helper class
			spanTree treeOut;
			treeOut.edge = prev;
			treeOut.cost = accumulate(dist.begin(),dist.end(),0);
			return treeOut;
			//If we have not reached all vertices, return a empty tree with cost -1
		}else if (Q.isEmpty()) return spanTree();
		
	}




int main(){
	
		//Testing read file constructor	
		std::ifstream inFile("cplusplus4c_homeworks_Homework3_SampleTestData_mst_data.txt"); 
		Graph newGraph(&inFile);
		std::cout<<newGraph<<std::endl;
		
		std::cout<<"Trial 1: "<<std::endl;
		spanTree result1 = newGraph.primMST(0); //Start growing the tree from node 0 by default;
		std::cout<<result1<<std::endl;
		
		std::cout<<"Trial 2: "<<std::endl;
		spanTree result2 = newGraph.primMST(19); //Also try start growing the tree from node 10;
		std::cout<<result2<<std::endl;
		

return 0;
}