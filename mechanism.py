"""
Author: Kevin J. Hwang
"""
import io
import math
import itertools
from preference import Preference
from profile import Profile

class Mechanism():
    """
    The parent class for all mechanisms. This class should not be constructed directly. All child
    classes are expected to contain the following variable(s).

    :ivar bool maximizeCandScore: True if the mechanism requires winners to maximize their score
        and if False otherwise.
    """

    def getWinners(self, profile):
        """
        Returns a list of all winning candidates given an election profile. This function assumes
        that getCandScoresMap(profile) is implemented for the child Mechanism class.
        
        :ivar Profile profile: A Profile object that represents an election profile.
        """
        
        candScores = self.getCandScoresMap(profile) 

        # Check whether the winning candidate is the candidate that maximizes the score or 
        # minimizes it.
        if self.maximizeCandScore == True:
            bestScore = max(candScores.values())
        else:
            bestScore = min(candScores.values())
        
        # Create a list of all candidates with the winning score and return it.
        winners = []
        for cand in candScores.keys():
            if candScores[cand] == bestScore:
                winners.append(cand)
        return winners

    def getRanking(self, profile):
        """
        Returns a list of lists that orders all candidates in tiers from best to worst given an 
        election profile. This function assumes that getCandScoresMap(profile) is implemented for 
        the child Mechanism class.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        # We generate a map that associates each score with the candidates that have that acore.
        candScoresMap = self.getCandScoresMap(profile) 
        reverseCandScoresMap = dict()
        for key, value in candScoresMap.items():
            if value not in reverseCandScoresMap.keys():
                reverseCandScoresMap[value] = [key]
            else:   
                reverseCandScoresMap[value].append(key)
        
        # We sort the scores by either decreasing order or increasing order.
        if self.maximizeCandScore == True:
            sortedCandScores = sorted(reverseCandScoresMap.keys(), reverse=True)
        else:
            sortedCandScores = sorted(reverseCandScoresMap.keys())
        
        # We put the candidates into our ranking based on the order in which their score appears
        ranking = []
        for candScore in sortedCandScores:
            currRanking = []
            for cand in reverseCandScoresMap[candScore]:
                currRanking.append(cand)
            ranking.append(currRanking)
        
        # Right now we return a list that contains the ranking list. This is for future extensions.
        results = []
        results.append(ranking)

        return results

class MechanismPosScoring(Mechanism):
    """
    The positional scoring mechanism. This class is the parent class for several mechanisms. This 
    can also be constructed directly. All child classes are expected to implement the
    getScoringVector() method.

    :ivar list<int> scoringVector: A list of integers (or floats) that give the scores assigned to
        each position in a ranking from first to last.
    """

    def __init__(self, scoringVector):
        self.maximizeCandScore = True
        self.scoringVector = scoringVector

    def isProfileValid(self, profile):
        elecType = profile.getElecType()
        if elecType != "soc" and elecType != "toc":
            return False
        return True

    def getScoringVector(self, profile):
        """
        Returns the scoring vector. This function is called by getCandScoresMap().

        :ivar Profile profile: A Profile object that represents an election profile.
        """
        
        # Check to make sure that the scoring vector contains a score for every possible rank in a
        # ranking.
        if len(self.scoringVector) != profile.numCands:
            print("ERROR: scoring vector is not the correct length")
            exit()
        
        return self.scoringVector

    def getCandScoresMap(self, profile):
        """
        Returns a dictonary that associates the integer representation of each candidate with the 
        score they recieved in the profile.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        # Currently, we expect the profile to contain complete ordering over candidates.
        elecType = profile.getElecType()
        if elecType != "soc" and elecType != "toc":
            print("ERROR: unsupported election type")
            exit()

        # Initialize our dictionary so that all candidates have a score of zero.
        candScoresMap = dict()
        for cand in profile.candMap.keys():
            candScoresMap[cand] = 0.0

        rankMaps = profile.getRankMaps()
        rankMapCounts = profile.getPreferenceCounts()
        scoringVector = self.getScoringVector(profile)

        # Go through the rankMaps of the profile and increment each candidates score appropriately.
        for i in range(0, len(rankMaps)):
            rankMap = rankMaps[i]
            rankMapCount = rankMapCounts[i]
            for cand in rankMap.keys():
                candScoresMap[cand] += scoringVector[rankMap[cand]-1]*rankMapCount
        
        return candScoresMap

    def getMov(self, profile):
        """
        Returns an integer that is equal to the margin of victory of the election profile, that is,
        the number of votes needed to be changed to change to outcome into a draw.

        :ivar Profile profile: A Profile object that represents an election profile.
        """
        import mov
        return mov.movPosScoring(profile, self.getScoringVector(profile))

class MechanismPlurality(MechanismPosScoring):
    """
    The plurality mechanism. This inherits from the positional scoring mechanism.
    """

    def __init__(self):
        self.maximizeCandScore = True

    def getScoringVector(self, profile):
        """
        Returns the scoring vector [1,0,0,...,0]. This function is called by getCandScoresMap() 
        which is implemented in the parent class.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        scoringVector = []
        scoringVector.append(1)
        for i in range(1, profile.numCands):
            scoringVector.append(0)
        return scoringVector

class MechanismVeto(MechanismPosScoring):
    """
    The veto mechanism. This inherits from the positional scoring mechanism.
    """

    def __init__(self):
        self.maximizeCandScore = True

    def getScoringVector(self, profile):
        """
        Returns the scoring vector [1,1,1,...,0]. This function is called by getCandScoresMap() 
        which is implemented in the parent class.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        numTiers = len(set(profile.getRankMaps()[0].values()))
        scoringVector = []
        for i in range(0, numTiers - 1):
            scoringVector.append(1)
        for i in range(numTiers - 1, profile.numCands):
            scoringVector.append(0)
        return scoringVector

class MechanismBorda(MechanismPosScoring):
    """
    The Borda mechanism. This inherits from the positional scoring mechanism.
    """

    def __init__(self):
        self.maximizeCandScore = True

    def getScoringVector(self, profile):
        """
        Returns the scoring vector [m-1,m-2,m-3,...,0] where m is the number of candidates in the 
        election profile. This function is called by getCandScoresMap() which is implemented in the
        parent class.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        scoringVector = []
        score = profile.numCands-1
        for i in range(0, profile.numCands):
            scoringVector.append(score)
            score -= 1
        return scoringVector

class MechanismKApproval(MechanismPosScoring):
    """
    The top-k mechanism. This inherits from the positional scoring mechanism.
    
    :ivar int k: The number of positions that recieve a score of 1. 
    """

    def __init__(self, k):
        self.maximizeCandScore = True
        self.k = k
    
    def getScoringVector(self, profile):
        """
        Returns a scoring vector such that the first k candidates recieve 1 point and all others 
        recive 0  This function is called by getCandScoresMap() which is implemented in the parent
        class.

        :ivar Profile profile: A Profile object that represents an election profile.
        """
        
        if self.k > profile.numCands:
            self.k = profile.numCands

        scoringVector = []
        for i in range(0, self.k):
            scoringVector.append(1)
        for i in range(self.k, profile.numCands):
            scoringVector.append(0)
        return scoringVector

class MechanismSimplifiedBucklin(Mechanism):
    """
    The simplified Bucklin mechanism.
    """

    def __init__(self):
        self.maximizeCandScore = False

    def getCandScoresMap(self, profile):
        """
        Returns a dictionary that associates integer representations of each candidate with their 
        Bucklin score.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        # Currently, we expect the profile to contain complete ordering over candidates.
        elecType = profile.getElecType()
        if elecType != "soc" and elecType != "toc":
            print("ERROR: unsupported profile type")
            exit()
        
        bucklinScores = dict()
        rankMaps = profile.getRankMaps()
        preferenceCounts = profile.getPreferenceCounts()
        for cand in profile.candMap.keys():

            # We keep track of the number of times a candidate is ranked in the first t positions.
            numTimesRanked = 0

            # We increase t in increments of 1 until we find t such that the candidate is ranked in the
            # first t positions in at least half the votes.
            for t in range(1, profile.numCands+1):
                for i in range(0, len(rankMaps)):        
                    if (rankMaps[i][cand] == t):
                        numTimesRanked += preferenceCounts[i]
                if numTimesRanked >= math.ceil(float(profile.numVoters)/2):
                    bucklinScores[cand] = t
                    break

        return bucklinScores

    def getMov(self, profile):
        """
        Returns an integer that is equal to the margin of victory of the election profile, that is,
        the number of votes needed to be changed to change to outcome into a draw.

        :ivar Profile profile: A Profile object that represents an election profile.
        """
        from . import mov
        return mov.movSimplifiedBucklin(profile)

class MechanismCopeland(Mechanism):
    """
    The Copeland mechanism.
    """

    def __init__(self, alpha):
        self.maximizeCandScore = True
        self.alpha = 0.5

    def getCandScoresMap(self, profile):
        """
        Returns a dictionary that associates integer representations of each candidate with their 
        Copeland score.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        # Currently, we expect the profile to contain complete ordering over candidates. Ties are
        # allowed however.
        elecType = profile.getElecType()
        if elecType != "soc" and elecType != "toc":
            print("ERROR: unsupported election type")
            exit()

        # Initialize each Copeland score as 0.0.
        copelandScores = dict()
        for cand in profile.candMap.keys():
            copelandScores[cand] = 0.0

        preferenceCounts = profile.getPreferenceCounts()

        # For each pair of candidates, calculate the number of votes in which one beat the other.
        wmgMap = profile.getWmg()
        for cand1, cand2 in itertools.combinations(wmgMap.keys(), 2):
            if cand2 in wmgMap[cand1].keys():
                if wmgMap[cand1][cand2] > 0:
                    copelandScores[cand1] += 1.0
                elif wmgMap[cand1][cand2] < 0:
                    copelandScores[cand2] += 1.0
            
                #If a pair of candidates is tied, we add alpha to their score for each vote.
                else:
                    copelandScores[cand1] += self.alpha
                    copelandScores[cand2] += self.alpha

        return copelandScores

class MechanismMaximin(Mechanism):
    """
    The maximin mechanism.
    """

    def __init__(self):
        self.maximizeCandScore = True

    def getCandScoresMap(self, profile):
        """
        Returns a dictionary that associates integer representations of each candidate with their 
        maximin score.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        # Currently, we expect the profile to contain complete ordering over candidates. Ties are
        # allowed however.
        elecType = profile.getElecType()
        if elecType != "soc" and elecType != "toc":
            print("ERROR: unsupported election type")
            exit()

        wmg = profile.getWmg()

        # Initialize the maximin score for each candidate as infinity.
        maximinScores = dict()
        for cand in wmg.keys():
            maximinScores[cand] = float("inf")
 
        # For each pair of candidates, calculate the number of times each beats the other.
        for cand1, cand2 in itertools.combinations(wmg.keys(), 2):
            if cand2 in wmg[cand1].keys():
                maximinScores[cand1] = min(maximinScores[cand1], wmg[cand1][cand2])
                maximinScores[cand2] = min(maximinScores[cand2], wmg[cand2][cand1])

        return maximinScores

class MechanismSchulze(Mechanism):
    """
    The Schulze mechanis, using the Floyd-Warshall algorithm.
    """

    def __init__(self):
        self.maximizeCandScore = True

    def computeStrongestPaths(self, profile, pairwisePreferences):
        """
        Returns a two-dimensional dictionary that associates every pair of candidates, cand1 and 
        cand2, with the strongest path from cand1 to cand2.

        :ivar Profile profile: A Profile object that represents an election profile.
        :ivar dict<int,dict<int,int>> pairwisePreferences: A two-dimensional dictionary that
            associates every pair of candidates, cand1 and cand2, with number of voters who prefer
            cand1 to cand2.
        """
        cands = profile.candMap.keys()
        numCands = len(cands)

        # Initialize the two-dimensional dictionary that will hold our strongest paths.
        strongestPaths = dict()
        for cand in cands:
            strongestPaths[cand] = dict()

        for i in range(1, numCands+1):
            for j in range(1, numCands+1):
                if (i == j):
                    continue
                if pairwisePreferences[i][j] > pairwisePreferences[j][i]:
                    strongestPaths[i][j] = pairwisePreferences[i][j]
                else:
                    strongestPaths[i][j] = 0

        for i in range(1, numCands+1):
            for j in range(1, numCands+1):
                if (i == j):
                    continue
                for k in range(1, numCands+1):
                    if (i == k or j == k):
                        continue
                    strongestPaths[j][k] = max(strongestPaths[j][k], min(strongestPaths[j][i], strongestPaths[i][k]))

        return strongestPaths

    def computePairwisePreferences(self, profile):
        """
        Returns a two-dimensional dictionary that associates every pair of candidates, cand1 and 
        cand2, with number of voters who prefer cand1 to cand2.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        cands = profile.candMap.keys()

        # Initialize the two-dimensional dictionary that will hold our pairwise preferences.
        pairwisePreferences = dict()
        for cand in cands:
            pairwisePreferences[cand] = dict()
        for cand1 in cands:    
            for cand2 in cands:
                if cand1 != cand2:
                    pairwisePreferences[cand1][cand2] = 0

        for preference in profile.preferences:
            wmgMap = preference.wmgMap
            for cand1, cand2 in itertools.combinations(cands, 2):
                
                # If either candidate was unranked, we assume that they are lower ranked than all
                # ranked candidates.
                if cand1 not in wmgMap.keys():
                    if cand2 in wmgMap.keys():
                        pairwisePreferences[cand2][cand1] += 1 * preference.count
                elif cand2 not in wmgMap.keys():
                    if cand1 in wmgMap.keys():
                        pairwisePreferences[cand1][cand2] += 1 * preference.count

                elif wmgMap[cand1][cand2] == 1:
                    pairwisePreferences[cand1][cand2] += 1 * preference.count
                elif wmgMap[cand1][cand2] == -1:
                    pairwisePreferences[cand2][cand1] += 1 * preference.count

        return pairwisePreferences

    def getCandScoresMap(self, profile):
        """
        Returns a dictionary that associates integer representations of each candidate with the
        number of other candidates for which her strongest path to the other candidate is greater
        than the other candidate's stronget path to her.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        cands = profile.candMap.keys()
        pairwisePreferences = self.computePairwisePreferences(profile)
        strongestPaths = self.computeStrongestPaths(profile, pairwisePreferences)

        # For each unranked candidate, find candidates where p[E,X] >= p[X,E], where X is all other
        # unranked candidates, using a variant of the Floyd-Warshall algorithm.
        ranking = dict()
        for cand in cands:
            ranking[cand] = 0

        winners = set()
        count = 0
        while len(winners) != len(cands):
            # store all winners found this round, as they will be tied in the ranking
            winners_temp = []
            for cand1 in cands:
                if cand1 in winners:
                    continue
                winner = True
                for cand2 in cands:
                    if cand2 in winners or cand1 == cand2:
                        continue
                    # if any unranked candidate beats this one, its not a winner
                    if strongestPaths[cand1][cand2] < strongestPaths[cand2][cand1]:
                        winner = False
                        break
                if winner:
                    winners_temp.append(cand1)
                    ranking[cand1] = len(cands) - count
            winners.update(winners_temp)
            count += 1
        return ranking

class MechanismSchulzeK(Mechanism):
    """
    The Schulze mechanism, using Kosaraju's algorithm.
    """

    def __init__(self):
        self.maximizeCandScore = True

    def visitDFS1(self, relation, search_order, visited, visit_ix, finish_order, finish_count):
        """
        Performs a depth-first search from a starting node and saves the pre-order of nodes visited
        Returns a count of how many nodes have been visited on the tree in total, and finish_order
        is updated to have the pre-order of nodes

        For use in the DFS of the original graph, in order to find the order to perform the second DFS in

        :ivar dict<int,dict<int,boolean>> relation: A two-dimensional dictionary that stores, for
            each candidate, whether they are pairwise preferred to another candidate or not
        :ivar list<int> search_order: A list that stores the order in which to go to new nodes in DFS
        :ivar dictionary<int,boolean> visited: A dictionary that stores whether each candidate's node
            has been visited or not
        :ivar int visit_ix: An int that represents which candidate node to visit
        :ivar list<int> finish_order: A list that stores the pre-order of nodes visited
        :ivar int finish_count: An integer that stores how many nodes have been visited in total
        """

        visited[visit_ix-1] = True

        for search_ix in range(len(search_order)):
            # new node to visit is determined by search_order
            probe_ix = search_order[search_ix]
            if relation[visit_ix][probe_ix] == True and visited[probe_ix-1] == False:
                # recursively continue DFS from the new node and update finish_count
                finish_count = self.visitDFS1(relation, search_order, visited, probe_ix, finish_order, finish_count)
        # store the finish order
        finish_order[finish_count] = visit_ix
        finish_count += 1
        return finish_count

    def visitDFS2(self, relation, search_order, visited, visit_ix, SCC):
        """
        Performs a depth-first search from a starting node and finds all nodes in the strongly connected
        component that the node is in
        Stores the discovered strongly connected component in SCC

        For use in the DFS of the transpose of the original graph, to find all SCCs

        :ivar dict<int,dict<int,boolean>> relation: A two-dimensional dictionary that stores, for
            each candidate, whether they are pairwise preferred to another candidate or not
        :ivar list<int> search_order: A list that stores the order in which to go to new nodes in DFS
        :ivar dictionary<int,boolean> visited: A dictionary that stores whether each candidate's node
            has been visited or not
        :ivar int visit_ix: An int that represents which candidate node to visit
        :ivar set<int> SCC: A set that stores all nodes in the strongly connected component
        """

        visited[visit_ix-1] = True
        # any node reached during DFS is in the same strongly connected component
        SCC.add(visit_ix)
        for search_ix in range(len(search_order)):
            # new node to visit is determined by search_order
            probe_ix = search_order[search_ix]
            if relation[visit_ix][probe_ix] == True and visited[probe_ix-1] == False:
                 # recursively continue DFS from the new node
                self.visitDFS2(relation, search_order, visited, probe_ix, SCC)


    #begins the DFS search
    #if second_DFS is true, DFS is done with visitDFS2 and returns tree_connects_any_prior_tree and tree
    #otherwise DFS is done with visitDFS1 and returns finish_order
    def startDFS(self, relation, search_order, N, second_DFS):
        """
        Initializes values for depth-first searches, and starts up appropriate DFS function
        Returns the per-order of the DFS if used for the first DFS, or the list of strongly connected
        components if used for the second DFS

        :ivar dict<int,dict<int,boolean>> relation: A two-dimensional dictionary that stores, for
            each candidate, whether they are pairwise preferred to another candidate or not
        :ivar list<int> search_order: A list that stores the order in which to go to new nodes in DFS
        :ivar int N: An int that represents the number of candidates
        :ivar boolean second_DFS: A boolean that is True if the second DFS function should be used,
            False otherwise
        """

        # intialize the visited list to False for all nodes
        visited = [False] * N

        # only care about finish order in first DFS
        if not second_DFS:
            finish_order = [-1] * N
            finish_count = 0        
        
        # only finding strongly connected components in second DFS
        if second_DFS:
            all_SCC = []

        # perform the appropriate DFS function
        for search_ix in range(N):
            root_ix = search_order[search_ix]
            if visited[root_ix-1] == False:
                if not second_DFS:
                    # with first DFS, generate pre-ordering of nodes
                    finish_count = self.visitDFS1(relation, search_order, visited, root_ix, finish_order, finish_count)
                else:
                    # with second DFS, find all strongly connected components
                    SCC = set()
                    self.visitDFS2(relation, search_order, visited, root_ix, SCC)
                    all_SCC.append(SCC)
        
        # return the appropriate value depending on which DFS was performed
        if not second_DFS:
            return finish_order
        else:
            return all_SCC

    def getCandScoresMap(self, profile):
        """
        Returns a dictionary that associates integer representations of each candidate with the
        number of other candidates for which her strongest path to the other candidate is greater
        than the other candidate's stronget path to her.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        cands = profile.candMap.keys()

        # set the number of candidates that need to be ranked by Kosaraju's algorithm
        # if set to len(cands), the full ranking will always be generated, otherwise
        # the algorithm will stop once at least this many candidates are ranked, but
        # may exceed this number
        top_rank = len(cands)

        # make sure that top_rank is less than 1 or greater than the number of candidates
        if top_rank < 1:
            top_rank = 1
        elif top_rank > len(cands):
            top_rank = len(cands)

        pairwise_preferences = self.computePairwisePreferences(profile)

        scores_map_prev = {}
        scores_map_final = {}

        while True:
            scores_map = self.kosaraju(pairwise_preferences,cands, [])

            # the first time the algorithm is run, copy the results into the final scores map
            if not scores_map_prev:
                scores_map_final = scores_map

            # if running Kosaraju's algorithm again didn't change the ranking, that means that
            # the candidates are actually fully ranked
            # this is needed to discriminate between unranked candidates and simply having
            # multiple candidates tied for last place
            if scores_map == scores_map_prev:
                break

            # add in the results of the scores_map
            for cand in cands:
                scores_map_final[cand] += scores_map[cand]

            # find the lowest score found so far, which is also the score of all unranked candidates
            bottom = min(set(scores_map.values()))

            # find all unranked candidates
            low_cands = set()
            for cand in cands:
                if scores_map[cand] == bottom:
                    low_cands.add(cand)
            # if only candidate is in the lowest tier, it is actually ranked by virtue of being last
            if len(low_cands) == 1:
                in_top = len(cands)
            # otherwise, the number of ranked candidates is the total number of candidates minus the unranked ones
            else:
                in_top = len(cands) - len(low_cands)

            # stop if enough candidates have been ranked
            if in_top >= top_rank - 1:
                break

            # delete all edges for already ranked candidates
            for cand1 in cands:
                if cand1 not in low_cands:
                    for cand2 in cands:
                        pairwise_preferences[cand1][cand2] = 0
                        pairwise_preferences[cand2][cand1] = 0

            # store the scores map for comparison in next iteration
            scores_map_prev = scores_map

        return scores_map_final
    
    def rankSchwartz(self, cands, pairwisePreferences):
        """
        Returns a dictionary that associates integer representations of each candidate with the
        number of other candidates for which her strongest path to the other candidate is greater
        than the other candidate's stronget path to her.

        Modified from the Floyd-Warshall implementation of Schulze, using parts of the
        computeStrongestPaths() and getCandScoresMap() functions
        For use in ranking candidates in the Schwartz set

        :ivar list<int> cands: A list of the candidates to rank
        :ivar dict<int,dict<int,int>> pairwisePreferences: A two-dimensional dictionary that
            associates every pair of candidates, cand1 and cand2, with number of voters who prefer
            cand1 to cand2.
        """
        numCands = len(cands)

        # Initialize the two-dimensional dictionary that will hold our strongest paths.
        strongestPaths = dict()
        for cand in cands:
            strongestPaths[cand] = dict()

        for i in cands:
            for j in cands:
                if (i == j):
                    continue
                if pairwisePreferences[i][j] > pairwisePreferences[j][i]:
                    strongestPaths[i][j] = pairwisePreferences[i][j]
                else:
                    strongestPaths[i][j] = 0

        for i in cands:
            for j in cands:
                if (i == j):
                    continue
                for k in cands:
                    if (i == k or j == k):
                        continue
                    strongestPaths[j][k] = max(strongestPaths[j][k], min(strongestPaths[j][i], strongestPaths[i][k]))

        # For each unranked candidate, find candidates where p[E,X] >= p[X,E], where X is all other
        # unranked candidates, using a variant of the Floyd-Warshall algorithm.
        ranking = dict()
        for cand in cands:
            ranking[cand] = 0

        winners = set()
        count = 0
        while len(winners) != len(cands):
            # store all winners found this round, as they will be tied in the ranking
            winners_temp = []
            for cand1 in cands:
                if cand1 in winners:
                    continue
                winner = True
                for cand2 in cands:
                    if cand2 in winners or cand1 == cand2:
                        continue
                    # if any unranked candidate beats this one, its not a winner
                    if strongestPaths[cand1][cand2] < strongestPaths[cand2][cand1]:
                        winner = False
                        break
                if winner:
                    winners_temp.append(cand1)
                    ranking[cand1] = len(cands) - count
            winners.update(winners_temp)
            count += 1
        return ranking

    def kosaraju(self, pairwise_preferences, cands, prev_maximal):
        """
        Returns a dictionary that associates integer representations of each candidate with the
        number of other candidates for which her strongest path to the other candidate is greater
        than the other candidate's stronget path to her.

        :ivar dict<int,dict<int,boolean>> pairwise_preferences: A two-dimensional dictionary that
            associates every pair of candidates, cand1 and cand2, with number of voters who prefer
            cand1 to cand2.
        :ivar list<int> cands: A list of the candidates to rank
        :ivar list<int> prev_maximal: A list of the candidates in the previously calculated Schwartz set
        """

        N = len(cands)

        # generate mappings of candidates that win pairwise comparisons, in the original graph
        # as well as the transpose
        relation = {}
        relation_transpose = {}
        for cand1 in cands:
            relation[cand1] = {}
            relation_transpose[cand1] = {}
            for cand2 in cands:
                cand1_beats_cand2 = False
                cand2_beats_cand1 = False
                if cand1 != cand2:
                    if pairwise_preferences[cand1][cand2] > pairwise_preferences[cand2][cand1]:
                        cand1_beats_cand2 = True
                    if pairwise_preferences[cand2][cand1] > pairwise_preferences[cand1][cand2]:
                        cand2_beats_cand1 = True
                relation[cand1][cand2] = cand1_beats_cand2
                relation_transpose[cand1][cand2] = cand2_beats_cand1

        #for first DFS, just search in order of candidate number
        search_order = range(1,len(cands)+1)

        #store the finish order from DFS
        finish_order = self.startDFS(relation, search_order, N, False)

        #for second DFS, the search order is the reverse of the finish order of first DFS
        finish_order.reverse()

        #find strongly connected components with second DFS
        all_SCC = self.startDFS(relation_transpose, finish_order, N, True)

        #find the candidates that are in the Schwartz set
        maximal = {}
        maximal_list = []

        for cand in cands:
            maximal[cand] = 0

        # for each SCC, determine if there are any incoming edges
        for SCC in all_SCC:
            is_maximal = True
            for cand1 in SCC:
                if not is_maximal:
                    break
                for cand2 in cands:
                    # if any candidate not in the SCC beats any candidate in the SCC,
                    # the SCC is not maximal
                    if relation[cand2][cand1] and not cand2 in SCC:
                        is_maximal = False
                        break
            if is_maximal:
                maximal_list.append(SCC)
                for cand in SCC:
                    maximal[cand] = 1

        # flatten the maximal list to produce a single list of all candidate in the Schwartz set
        all_maximal = [val for sublist in maximal_list for val in sublist]

        # check whether member of the Schartz set beat any other members of the Schwartz set
        no_defeats = True
        for cand1 in all_maximal:
            for cand2 in all_maximal:
                if relation[cand2][cand1]:
                    no_defeats = False
                    break

        # if there were no defeats in the Schwartz set, the winners have been found
        if no_defeats:
            # if there is a previously calculated Schwartz set, use the Floyd-Warshall algorithm
            # to produce a ranking of the Schwartz set
            if prev_maximal:
                ranking = self.rankSchwartz(prev_maximal, pairwise_preferences)
                for cand in prev_maximal:
                    maximal[cand] += ranking[cand]
            return maximal

        # if there were defeats in the Schwartz set, perform sequential dropping of the weakest defeats
        weakest_defeat = -1

        # first find the value of the weakest defeat
        for cand1 in all_maximal:
            for cand2 in all_maximal:
                if relation[cand1][cand2] and (weakest_defeat > pairwise_preferences[cand1][cand2] or weakest_defeat == -1):
                    weakest_defeat = pairwise_preferences[cand1][cand2]

        # if the weakest defeat is not a tie, drop the weakest defeat
        if weakest_defeat >= 0:
            for cand1 in all_maximal:
                for cand2 in all_maximal:
                    # set all defeats equal to the weakest defeat to be ties
                    if relation[cand1][cand2] and pairwise_preferences[cand1][cand2] == weakest_defeat:
                        pairwise_preferences[cand1][cand2] = 0
                        pairwise_preferences[cand2][cand1] = 0

            # recursively calculate the next Schwartz set
            next_maximal = self.kosaraju(pairwise_preferences,cands,all_maximal)
            # add the results of the new Schwartz set to what has already been calculated
            for cand in next_maximal:
                maximal[cand] += next_maximal[cand]

        return maximal

    def computePairwisePreferences(self, profile):
        """
        Returns a two-dimensional dictionary that associates every pair of candidates, cand1 and 
        cand2, with number of voters who prefer cand1 to cand2.

        :ivar Profile profile: A Profile object that represents an election profile.
        """

        cands = profile.candMap.keys()

        # Initialize the two-dimensional dictionary that will hold our pairwise preferences.
        pairwisePreferences = dict()
        for cand in cands:
            pairwisePreferences[cand] = dict()
        for cand1 in cands:    
            for cand2 in cands:
                if cand1 != cand2:
                    pairwisePreferences[cand1][cand2] = 0

        for preference in profile.preferences:
            wmgMap = preference.wmgMap
            for cand1, cand2 in itertools.combinations(cands, 2):
                
                # If either candidate was unranked, we assume that they are lower ranked than all
                # ranked candidates.
                if cand1 not in wmgMap.keys():
                    if cand2 in wmgMap.keys():
                        pairwisePreferences[cand2][cand1] += 1 * preference.count
                elif cand2 not in wmgMap.keys():
                    if cand1 in wmgMap.keys():
                        pairwisePreferences[cand1][cand2] += 1 * preference.count

                elif wmgMap[cand1][cand2] == 1:
                    pairwisePreferences[cand1][cand2] += 1 * preference.count
                elif wmgMap[cand1][cand2] == -1:
                    pairwisePreferences[cand2][cand1] += 1 * preference.count

        return pairwisePreferences


def getKendallTauScore(myResponse, otherResponse):
    """
    Returns the Kendall Tau Score
    """
    # variables
    kt = 0
    list1 = myResponse.values()
    list2 = otherResponse.values()
    
    if len(list1) <= 1:
        return kt

    #runs through list1
    for itr1 in range(0, len(list1) - 1):
        #runs through list2
        for itr2 in range(itr1 + 1, len(list2)):
            # checks if there is a discrepancy. If so, adds
            if ((list1[itr1] > list1[itr2]
                and list2[itr1] < list2[itr2])
            or (list1[itr1] < list1[itr2]
                and list2[itr1] > list2[itr2])):

                kt += 1
    # normalizes between 0 and 1
    kt = (kt * 2) / (len(list1) * (len(list1) - 1))

    #returns found value
    return kt
