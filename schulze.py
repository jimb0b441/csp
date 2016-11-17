#performs Kosaraju's algorithms to determine Schulze winners
#returns dictionary of candidates to to their Schulze ranking. All candidates not in Schwartz set are ranked as 0
def getCandScoresMap(self, profile):
    cands = profile.candMap.keys()
    N = len(cands)
    pairwise_preferences = self.computePairwisePreferences(profile)
    
    #generate mapping of candidates that win pairwise comparisons
    relation = {}
    for cand1 in cands:
        relation[cand1] = {}
        for cand2 in cands:
            visit = False
            if cand1 != cand2 and pairwise_preferences[cand1][cand2] > pairwise_preferences[cand2][cand1]:
                visit = True
            relation[cand1][cand2] = visit
    
    #for first DFS, just search in order of candidate number
    search_order = range(1,len(cands+1))
    #store the finish order from DFS
    finish_order = startDFS(relation, search_order, N, False)
    
    #for second DFS, the search order is the reverse of the finish order of first DFS
    search_order2 = finish_order.reverse()
    #new relation is the transpose of the first
    relation2 = {}
    for cand1 in cands:
        relation2[cand1] = {}
        for cand2 in cands:
            relation2[cand1][cand2] = relation[cand2][cand1]
    
    #find strongly connected components with second DFS
    tree_connects_any_prior_tree, tree = startDFS(relation2, search_order2, N, True)
    
    #find the candidates that are in the Schwartz set
    maximal = {}
    maximal_list = []
    for cand in cands:
        if tree_connects_any_prior_tree[tree[cand]] == False:
            maximal[cand] = 1
            maximal_list.append(cand)
        else:
            maximal[cand] = 0
    
    #remove ambiguity
    while len(maximal_list) > 1:
        weakest_defeat = N+1
        strongest_defeat = -1
        for cand1 in maximal_list:
            for cand2 in maximal_list:
                if cand1 != cand2:
                    defeat = pairwise_preferences[cand1][cand2]
                    if defeat > pairwise_preferences[cand2][cand1]:
                        if defeat < weakest_defeat:
                            weakest_defeat = defeat
                        if defeat > strongest_defeat:
                            strongest_defeat = defeat
        if weakest_defeat == strongest_defeat:
            break
        new_maximal = maximal_list
        for cand1 in maximal_list:
            for cand2 in maximal_list:
                if cand1 != cand2:
                    defeat = pair_wise_preferences[cand1][cand2]
                    if defeat > pairwise_preferences[cand2][cand1] and defeat == weakest_defeat:
                        new_maximal.remove(cand1)
        maximal_list = new_maximal
        for cand in maximal_list:
            maximal[cand] += 1
        
    return maximal
    
#begins the DFS search
#if second_DFS is true, DFS is done with visitDFS2 and returns tree_connects_any_prior_tree and tree
#otherwise DFS is done with visitDFS1 and returns finish_order
def startDFS(relation, search_order, N, second_DFS):
    visited = [False] * N
    tree = [0] * N
    tree_count = 0

    #only care about finish order in first DFS
    if not second_DFS:
        finish_order = [-1] * N
        finish_count = 0        
    
    #only finding strongly connected components in second DFS
    if second_DFS:
        tree_connects_any_prior_tree = [False] * N
    
    for search_ix in range(N):
        root_ix = search_order[search_ix]
        if visited[root_ix] == False:
            if not second_DFS:
                finish_count = visitDFS1(relation, visited, search_order, root_ix, finish_order, finish_count, tree, tree_count)
            else:
                visitDFS2(relation, visited, search_order, root_ix, tree, tree_connects_any_prior_tree, tree_count)
            tree_count += 1
            
    if not second_DFS:
        return finish_order
    else:
        return tree_connects_any_prior_tree, tree

#DFS for the first search
#returns finish_count
def visitDFS1(relation, search_order, visited, visit_ix, finish_order, finish_count, tree, tree_count):
    tree[visit_ix] = tree_count
    visited[visit_ix] = True
    
    for search_ix in range(N):
        probe_ix = search_order[search_ix]
        if relation[visit_ix][probe_ix] == True and visited[probe_ix] == False:
            finish_count = visitDFS1(relation, search_order, visited, probe_ix, finish_order, finish_count, tree, tree_count)

    finish_count += 1
    finish_order[finish_count] = visit_ix
    return finish_count

#DFS for second search
#may be able to simplify further, see note below pseudo-code
def visitDFS2(relation, search_order, visited, visit_ix, tree, tree_connects_any_prior_tree, tree_count):
    tree[visit_ix] = tree_count
    visited[visit_ix] = True
    
    for search_ix in range(N):
        probe_ix = search_order[search_ix]
        if relation[visit_ix][probe_ix] == True:
            if visited[probe_ix] == False:
                visitDFS2(relation, search_order, visited, probe_ix, tree, tree_connects_any_prior_tree, tree_count)
            else:
                if tree[prob_ix] < tree[tree_count]:
                    tree_connects_any_prior_tree[tree_count] = True
                    
    
