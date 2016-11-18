#to use for testing implementation of schulze

from profile import Profile
from preference import Preference
import mechanism 
import sys
import numpy as np
import time

def genPref(ranking, num_votes):
	wmgMap = dict()
	for i in range(1,len(ranking)+1):
		wmgMap[i] = dict()
	for i in ranking:
		for j in ranking:
			if i != j:
				if ranking.index(i) < ranking.index(j):
					wmgMap[i][j] = 1
				else:
					wmgMap[i][j] = -1
	return Preference(wmgMap, num_votes)

def Main():
	icand = int(sys.argv[1]) # no. of candidates
	iranks = int(sys.argv[2]) # no. of rankings

	preferences = []
	for i in range(iranks):
		preferences.append(genPref((np.random.permutation(icand)+1).tolist(), np.random.randint(iranks)+1))

	candMap = dict()
	for i in range(icand):
		candMap[i+1] = i+1

	profile = Profile(candMap, preferences)

	start = time.time()
	mechanism1 = mechanism.MechanismSchulzeFW()
	rankingFW = mechanism1.getRanking(profile)
	end = time.time()
	print()
	print('Floyd-Warshall elapsed time: ', end-start)

	start = time.time()
	mechanism2 = mechanism.MechanismSchulzeK()
	rankingK = mechanism2.getRanking(profile)
	end = time.time()
	print()
	print("Kosaraju elapsed time: ", end-start)

if __name__ == '__main__':
	Main()
