#include "Graph.h"

void Graph::Dfs(int v) 
{
	used[v] = true;
	comp.push_back(v);
	for (size_t i = 0; i < g[v].size(); ++i) 
	{
		int to = g[v][i];
		if (!used[to])
			Dfs(to);
	}
}

void Graph::FindComps() 
{
	for (int i = 0; i < n; ++i)
		used[i] = false;
	for (int i = 0; i < n; ++i)
		if (!used[i]) 
		{
			comp.clear();
			Dfs(i);

			cout << "Component:";
			for (size_t j = 0; j < comp.size(); ++j)
				cout << ' ' << comp[j];
			cout << endl;
		}
}

