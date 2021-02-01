#pragma once

#include "stdafx.h"
#include <stack>

class Graph
{
private:
	int n;
	vector<vector<int>> g;
	vector<bool> used;
	vector<int> comp;

public:
	void init();
	void Dfs(int v);
	void FindComps();
};
