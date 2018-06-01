
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
using namespace std;

class TopAdjList;
class AdjMatrix;
class EdgeList;

class DSU {
public:
	DSU(size_t size) {
		_dsu = vector<size_t>(size);
		for (size_t i = 0; i <_dsu.size(); ++i)
			_dsu[i] = i;
		_size = size;
	}
	size_t find(size_t v)
	{
		if (v == _dsu[v])
			return v;
		return _dsu[v] = find(_dsu[v]);
	}
	void unite(size_t first, size_t second)
	{
		first = find(first);
		second = find(second);
		if (rand() & 1)
			swap(first, second);
		if (first != second) {
			_dsu[first] = second;
			--_size;
		}
	}
	size_t size()
	{
		return _size;
	}

	vector<set<size_t>> sets()
	{
		map<size_t, set<size_t>> sets = map<size_t, set<size_t>>();
		for (size_t i = 0; i<_size; ++i)
			sets[i].insert(find(i));
		vector<set<size_t>> result = vector<set<size_t>>();
		for (auto it : sets)
			result.push_back(it.second);
		return result;
	}
private:
	vector<size_t> _dsu;
	size_t _size;
};

class GraphPresentation {
public:
	char pres = 'U';
	int N;
	bool isWeighted;
	bool isOriented;
	virtual void readGraph(istream &file) {
		file >> isOriented >> isWeighted;
	}
	virtual void writeGraph(ostream & os) { os << isOriented << ' ' << isWeighted << endl; }
	virtual void addEdge(int from, int to, int weight) = 0;
	virtual int changeEdge(int from, int to, int newWeight) = 0;
	virtual void removeEdge(int from, int to) = 0;
	virtual GraphPresentation * transformToAdjList() = 0;
	virtual GraphPresentation * transformToAdjMatrix() = 0;
	virtual GraphPresentation * transformToListOfEdges() = 0;
	virtual GraphPresentation * getSpaingTreePrima() = 0;
	virtual GraphPresentation * getSpaingTreeKruscal() = 0;
	virtual GraphPresentation * getSpaingTreeBoruvka() = 0;

};

class Graph {
public:
	int N;
	GraphPresentation * presentation;

	Graph() { N = 0; }
	Graph(int _N, GraphPresentation * pres) { N = _N; presentation = pres; }
	void readGraph(string fileName);
	void writeGraph(string fileName);
	void addEdge(int from, int to, int weight);
	void removeEdge(int from, int to);
	int changeEdge(int from, int to, int newWeight);
	void transformToAdjList();
	void transformToAdjMatrix();
	void transformToListOfEdges();
	Graph getSpaingTreePrima();
	Graph getSpaingTreeKruscal();
	Graph getSpaingTreeBoruvka();
};



class AdjMatrix : public GraphPresentation {
public:
	char pres = 'C';
	vector<vector<int>> matrix;	// представление этого класса

	AdjMatrix() { this->N = 0; }
	AdjMatrix(int _N) { this->N = _N; }
	AdjMatrix(int _N, bool _isOriented, bool _isWeighted, vector<vector<int>> adjMatrix) { this->N = _N; this->isOriented = _isOriented; this->isWeighted = _isWeighted; this->matrix = adjMatrix; }

	void readGraph(istream &file) {
		matrix.resize(N);
		for (int i = 0; i < N; i++) {
			matrix[i].resize(N);
			for (int j = 0; j < N; j++) {
				file >> matrix[i][j];
				/*if (!this->isOriented)
				this->_matrix[j][i] = weight;*/
			}
		}
	}

	void writeGraph(ostream & os) {
		os << this->pres << ' ';
		cout << this->pres << ' ';
		os << this->N << endl;
		GraphPresentation::writeGraph(os);
		for (int i = 0; i < this->N; i++)
			for (int j = 0; j < this->N; j++)
				if (j < N - 1)
					os << this->matrix[i][j] << ' ';
				else
					os << this->matrix[i][j] << endl;
	}

	void addEdge(int from, int to, int weight) {
		if (!this->isWeighted)
			weight = 1;
		matrix[from][to] = weight;
		if (!this->isOriented)
			matrix[to][from] = weight;
	}

	void removeEdge(int from, int to) {
		this->matrix[from][to] = 0;
		if (!this->isOriented)
			this->matrix[to][from] = 0;
	}

	int changeEdge(int from, int to, int newWeight) {
		if (!this->isWeighted)
			newWeight = 1;
		int oldWeight = this->matrix[from][to];
		this->matrix[from][to] = newWeight;
		if (!this->isOriented)
			this->matrix[to][from] = newWeight;
		return oldWeight;
	}

	GraphPresentation * transformToAdjList();

	GraphPresentation * transformToAdjMatrix() { return this; }

	GraphPresentation * transformToListOfEdges();

	GraphPresentation * getSpaingTreePrima() { return this; }

	GraphPresentation * getSpaingTreeKruscal() { return this; }

	GraphPresentation * getSpaingTreeBoruvka() { return this; }
};

class TopAdjList : public GraphPresentation {
public:
	char pres = 'L';
	vector<map<int, int>> topAdjList;	// представление этого класса

	TopAdjList(int _N) { this->N = _N; }
	TopAdjList(int _N, bool _isOriented, bool _isWeighted, vector<map<int, int>> adjList) { this->N = _N; this->isOriented = _isOriented; this->isWeighted = _isWeighted; this->topAdjList = adjList; }

	void readGraph(istream &file) {
		GraphPresentation::readGraph(file);
		int neighbor;
		int weight;
		char notused;	// пробел
		for (int i = 0; i < this->N; i++) {
			file >> neighbor;
			neighbor--;
			if (this->isWeighted)
				file >> notused >> weight;
			else
				weight = 1;
			if (this->isOriented)
				topAdjList[i][neighbor] = weight;
			else
			{
				topAdjList[i][neighbor] = weight;
				topAdjList[neighbor][i] = weight;
			}
		}

	}

	void writeGraph(ostream & os) {
		os << this->pres << ' ';
		//cout << this->pres << ' ';
		os << this->N << endl;
		//cout << this->N << endl;
		GraphPresentation::writeGraph(os);

		for (int i = 0; i < this->N; i++)
		{
			for (auto pair : topAdjList[i])
			{
				if (this->isWeighted) {
					os << pair.first + 1 << ' ' << pair.second << ' ';
					//cout << pair.first + 1 << ' ' << pair.second << ' ';
				}
				else {
					os << pair.first + 1 << ' ';
					//cout << pair.first + 1 << ' ';
				}
				os << endl;
				//cout << endl;
			}
		}
		//cout << endl;
	}

	void addEdge(int from, int to, int weight) {
		topAdjList[from].insert(pair<int, int>(to, weight));
		if (!this->isOriented)
			topAdjList[to].insert(pair<int, int>(from, weight));
	}

	void removeEdge(int from, int to) {
		topAdjList[from].erase(to);
		if (!this->isOriented)
			topAdjList[to].erase(from);
	}

	int changeEdge(int from, int to, int newWeight) {
		if (!this->isWeighted)
			newWeight = 1;

		int oldWeight = topAdjList[from][to];
		topAdjList[from][to] = newWeight;
		if (!this->isOriented)
			topAdjList[to][from] = newWeight;
		return oldWeight;
	}

	GraphPresentation * transformToAdjList() { return this; }

	GraphPresentation * transformToAdjMatrix() {
		if (this->pres == 'L') {
			auto& adjList = this->topAdjList;
			vector<vector<int>> adjMatrix = vector<vector<int>>(this->N, vector<int>(this->N));
			for (int k = 0; k < adjList.size(); k++)
				for (auto it : adjList[k])
					adjMatrix[k][it.first] = it.second;
			return new AdjMatrix(this->N, this->isOriented, this->isWeighted, adjMatrix);
		}
	}

	GraphPresentation * transformToListOfEdges();

	GraphPresentation * getSpaingTreePrima() {
		int vertices = this->N;
		int edges = 0;
		int next_vertex;
		vector<int> actualDistance = vector<int>(N, INT32_MAX);
		vector<int> parent = vector<int>(N, -1);
		vector<bool> visited = vector<bool>(N, false);

		auto comp = [](pair<int, int> a, pair<int, int> b) {
			if (a.second == b.second)
				return a.first < b.first;
			else
				return a.second < b.second;
		};
		auto possibleDistanceToTree = set<pair<int, int>, decltype(comp)>(comp);

		auto temp = vector<map<int, int>>(N);
		GraphPresentation* spanning_tree = new TopAdjList(N, this->isOriented, this->isWeighted, temp);

		for (size_t i = 1; i < vertices; i++)
			possibleDistanceToTree.insert(make_pair(i, INT32_MAX));

		next_vertex = 0;
		parent[0] = 0;
		visited[0] = true;
		while (edges < N - 1) {
			for (auto v : topAdjList[next_vertex])
				if ((parent[v.first] == -1 || actualDistance[v.first] > v.second) & !visited[v.first]) {
					actualDistance[v.first] = v.second;
					parent[v.first] = next_vertex;
					possibleDistanceToTree.insert(make_pair(v.first, v.second));
				}
			auto it = *possibleDistanceToTree.begin();
			possibleDistanceToTree.erase(possibleDistanceToTree.begin());
			while (actualDistance[it.first] != it.second) {
				it = *possibleDistanceToTree.begin();
				possibleDistanceToTree.erase(possibleDistanceToTree.begin());
			}
			next_vertex = it.first;
			spanning_tree->addEdge(parent[next_vertex], next_vertex, it.second);
			visited[next_vertex] = true;
			++edges;
		}
		return spanning_tree;
	}

	GraphPresentation * getSpaingTreeKruscal() { return this; }

	GraphPresentation * getSpaingTreeBoruvka() { return this; }

};

class EdgeList : public GraphPresentation {
public:
	char pres = 'E';
	map<pair<int, int>, int> edgeList;	// представление этого класса
	int M;

	EdgeList(int _N) { this->N = _N; }
	EdgeList(int _N, bool _isOriented, bool _isWeighted, int size, map<pair<int, int>, int> _edgeList) { this->N = _N; this->isOriented = _isOriented; this->isWeighted = _isWeighted; this->M = size; this->edgeList = _edgeList; }
	void readGraph(istream &file) {
		int _M;

		file >> _M;	this->M = _M;
		GraphPresentation::readGraph(file);

		int node1, node2, weight;
		for (int i = 0; i < _M; i++) {
			file >> node1 >> node2;
			node1--;
			node2--;
			if (this->isWeighted)
				file >> weight;
			else
				weight = 1;
			this->addEdge(node1, node2, weight);
		}

	}

	void writeGraph(ostream & os) {
		os << this->pres << ' ';
		//cout << this->pres << ' ';
		os << N << ' ' << M << endl;
		//cout << N << ' ' << M << endl;
		GraphPresentation::writeGraph(os);

		for (auto pair : edgeList)
		{
			auto nodes = pair.first;
			auto weight = pair.second;
			if (this->isWeighted) {
				os << nodes.first + 1 << ' ' << nodes.second + 1 << ' ' << weight << endl;
				//cout << nodes.first + 1 << ' ' << nodes.second + 1 << ' ' << weight << endl;
			}
			else {
				os << nodes.first + 1 << ' ' << nodes.second + 1 << endl;
				//cout << nodes.first + 1 << ' ' << nodes.second + 1 << endl;
			}
		}
		//cout << endl;
	}

	void addEdge(int from, int to, int weight) {
		edgeList[pair<int, int>(from, to)] = weight;
	}

	void removeEdge(int from, int to) {
		M--;
		edgeList.erase(pair<int, int>(from, to));
	}

	int changeEdge(int from, int to, int newWeight) {
		if (!this->isWeighted)
			newWeight = 1;

		int oldWeight = edgeList[make_pair(from, to)];
		edgeList[pair<int, int>(from, to)] = newWeight;
		return oldWeight;
	}

	GraphPresentation * transformToAdjList() {
		if (this->pres == 'E') {
			auto& _edgeList = this->edgeList;
			vector<map<int, int>> adjList = vector<map<int, int>>(this->N);
			for (auto it : _edgeList)
			{
				auto& weight = it.second;
				auto& vertices = it.first;
				adjList[vertices.first][vertices.second] = weight;
				if (!this->isOriented)
					adjList[vertices.second][vertices.first] = weight;
			}
			return new TopAdjList(this->N, this->isOriented, this->isWeighted, adjList);
		}
	}

	GraphPresentation * transformToAdjMatrix() {
		if (this->pres == 'E') {
			auto& edgeList = this->edgeList;
			vector<vector<int>> adjMatrix = vector<vector<int>>(this->N, vector<int>(this->N));
			for (auto it : edgeList)
			{
				auto& weight = it.second;
				auto& vertices = it.first;
				adjMatrix[vertices.first][vertices.second] = weight;
				if (!this->isOriented)
					adjMatrix[vertices.second][vertices.first] = weight;
			}
			return new AdjMatrix(this->N, this->isOriented, this->isWeighted, adjMatrix);
		}
	}

	GraphPresentation * transformToListOfEdges() { return this; }

	GraphPresentation * getSpaingTreeKruscal() {
		int vertices = N;
		int edges = 0;
		GraphPresentation* spanning_tree = new EdgeList(vertices, this->isOriented, this->isWeighted, edges, *(new map<pair<int, int>, int>()));
		DSU dsu = DSU(vertices);
		auto comp = [](pair<pair<int, int>, int> a, pair<pair<int, int>, int> b) {
			if (a.second == b.second)
				if (a.first.first == b.first.first)
					if (a.first.second == b.first.second)
						return false;
					else
						return a.first.second < b.first.second;
				else
					return a.first.first < b.first.first;
			else
				return a.second < b.second; };
		auto sortedEdges = set<pair<pair<int, int>, int>, decltype(comp)>(comp);
		for (auto e : this->edgeList)
			sortedEdges.insert(make_pair(make_pair(e.first.first, e.first.second), e.second));
		while (edges < vertices - 1)
		{
			cout << endl;
			auto next_edge = *sortedEdges.begin();
			sortedEdges.erase(sortedEdges.begin());
			size_t from = next_edge.first.first, to = next_edge.first.second, weight = next_edge.second;
			if (dsu.find(from) != dsu.find(to))
			{
				dsu.unite(from, to);
				spanning_tree->addEdge(from, to, weight);
				++edges;
			}
		}
		return spanning_tree;
	}

	GraphPresentation * getSpaingTreeBoruvka() {
		int vertices = this->N;
		int edges = 0;
		auto temp = map<pair<int, int>, int>();
		GraphPresentation* spanning_tree = new EdgeList(vertices, this->isOriented, this->isWeighted, edges, temp);
		DSU dsu = DSU(vertices);
		vector<bool> foundEdges = vector<bool>(vertices, false);
		auto minEdges = vector<pair<pair<int, int>, int>>(vertices);
		while (dsu.size() != 1)
		{
			std::fill(foundEdges.begin(), foundEdges.end(), false);
			for (auto e : this->edgeList)
			{
				size_t comp1 = dsu.find(e.first.first), comp2 = dsu.find(e.first.second);
				if (comp1 == comp2)
					continue;
				int weight = e.second;
				if (!foundEdges[comp1] || minEdges[comp1].second > weight)
				{
					minEdges[comp1] = e;
					foundEdges[comp1] = true;
				}
				if (!foundEdges[comp2] || minEdges[comp2].second > weight)
				{
					minEdges[comp2] = e;
					foundEdges[comp2] = true;
				}
			}
			for (auto e : minEdges)
			{
				int from = e.first.first, to = e.first.second;
				int weight = e.second;
				if (dsu.find(from) == dsu.find(to))
					continue;
				spanning_tree->addEdge(from, to, weight);
				++edges;
				dsu.unite(from, to);
			}
		}
		return spanning_tree;
	}

	GraphPresentation * getSpaingTreePrima() { return this; }
};

GraphPresentation * AdjMatrix::transformToListOfEdges() {
	if (this->pres == 'C') {
		auto& adjMatrix = this->matrix;
		map<pair<int, int>, int> edgeList;
		for (size_t k = 0; k < adjMatrix.size(); k++)
			for (size_t m = 0; m < adjMatrix[k].size(); m++)
				if (adjMatrix[k][m] != 0)
					if (this->isOriented || edgeList.find(make_pair(m, k)) == edgeList.end())
						edgeList[make_pair(k, m)] = adjMatrix[k][m];
		return new EdgeList(this->N, this->isOriented, this->isWeighted,
			edgeList.size(), edgeList);
	}
}

GraphPresentation * AdjMatrix::transformToAdjList() {
	if (this->pres == 'C') {
		auto& adjMatrix = this->matrix;
		vector<map<int, int>> adjList = vector<map<int, int>>(this->N);
		for (int k = 0; k < adjMatrix.size(); k++)
			for (int m = 0; m < adjMatrix[k].size(); m++)
				if (adjMatrix[k][m] != 0)
					adjList[k][m] = adjMatrix[k][m];
		return new TopAdjList(this->N, this->isOriented, this->isWeighted, adjList);
	}
}

GraphPresentation * TopAdjList::transformToListOfEdges() {
	if (this->pres == 'L') {
		auto& adjList = this->topAdjList;
		map<pair<int, int>, int> edgeList;
		for (size_t k = 0; k < adjList.size(); k++)
			for (auto it : adjList[k])
				if (this->isOriented || edgeList.find(make_pair(it.first, k)) == edgeList.end())
					edgeList[make_pair(k, it.first)] = it.second;
		return new EdgeList(this->N, this->isOriented, this->isWeighted, edgeList.size(), edgeList);
	}
}

void Graph::readGraph(string fileName)
{
	ifstream file(fileName);
	//file.open(fileName);
	int _N;
	char presentation;
	file >> presentation >> _N;

	this->presentation = new AdjMatrix();
	this->presentation->N = _N;
	this->N = _N;
	switch (presentation)
	{
	case 'C':
		this->presentation = new AdjMatrix(_N);
		break;
	case 'L':
		this->presentation = new TopAdjList(_N);
		break;
	case 'E':
		this->presentation = new EdgeList(_N);
		break;
	}

	this->presentation->readGraph(file);
	file.close();
}

void Graph::writeGraph(string fileName) {
	ofstream file;
	file.open(fileName);
	//file << this->presentation->pres << ' ';
	//cout << this->presentation->pres << ' ';
	this->presentation->writeGraph(file);
	file.close();
}

void Graph::addEdge(int from, int to, int weight) {
	this->presentation->addEdge(from, to, weight);
}

void Graph::removeEdge(int from, int to) {
	this->presentation->removeEdge(from, to);
}

int Graph::changeEdge(int from, int to, int newWeight) {
	return this->presentation->changeEdge(from, to, newWeight);
}

void Graph::transformToAdjList() {
	if (this->presentation->pres == 'L')
		return;
	GraphPresentation * newPresentation = this->presentation->transformToAdjList();
	this->presentation = newPresentation;
}

void Graph::transformToAdjMatrix() {
	if (this->presentation->pres == 'L')
		return;
	GraphPresentation * newPresentation = this->presentation->transformToAdjMatrix();
	this->presentation = newPresentation;
}

void Graph::transformToListOfEdges() {
	if (this->presentation->pres == 'L')
		return;
	GraphPresentation * newPresentation = this->presentation->transformToListOfEdges();
	this->presentation = newPresentation;

}

Graph Graph::getSpaingTreePrima() {

	this->transformToAdjList();
	GraphPresentation * graphPres = this->presentation->getSpaingTreePrima();
	auto result = Graph(this->N, graphPres);

	return result;
}

Graph Graph::getSpaingTreeKruscal() {

	this->transformToListOfEdges();
	GraphPresentation* graphPres = this->presentation->getSpaingTreeKruscal();
	auto result = Graph(this->N, graphPres);

	return result;
}

Graph Graph::getSpaingTreeBoruvka() {

	this->transformToListOfEdges();
	GraphPresentation* graphPres = this->presentation->getSpaingTreeBoruvka();
	auto result = Graph(this->N, graphPres);

	return result;
}



int main()
{
	Graph g;
	g.readGraph("input.txt"); g.writeGraph("output.txt");
	//Graph gg=g.getSpaingTreeBoruvka();
	Graph gg = g.getSpaingTreeKruscal(); gg.writeGraph("output.txt");
	// Graph gg=g.getSpaingTreePrima();
	gg.transformToAdjList(); 
	gg.writeGraph("output.txt");

	getchar();
	return 0;
}

