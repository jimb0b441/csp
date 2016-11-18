#to use for testing implementation of schulze

from profile import Profile
from preference import Preference
import mechanism 

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


preferences = []

#c>d>b>a, with c, d, and b in a cycle (schwartz set)
#tie breaking needed to rank the schwartz set
preferences.append(genPref([1,2,3,4], 0))
preferences.append(genPref([1,2,4,3], 0))
preferences.append(genPref([1,3,2,4], 0))
preferences.append(genPref([1,3,4,2], 0))
preferences.append(genPref([1,4,2,3], 0))
preferences.append(genPref([1,4,3,2], 0))
preferences.append(genPref([2,1,4,3], 0))
preferences.append(genPref([2,1,3,4], 0))
preferences.append(genPref([2,3,4,1], 2))
preferences.append(genPref([2,3,1,4], 0))
preferences.append(genPref([2,4,3,1], 0))
preferences.append(genPref([2,4,1,3], 0))
preferences.append(genPref([3,1,2,4], 0))
preferences.append(genPref([3,1,4,2], 0))
preferences.append(genPref([3,2,1,4], 0))
preferences.append(genPref([3,2,4,1], 0))
preferences.append(genPref([3,4,1,2], 0))
preferences.append(genPref([3,4,2,1], 4))
preferences.append(genPref([4,1,3,2], 0))
preferences.append(genPref([4,1,2,3], 0))
preferences.append(genPref([4,2,3,1], 3))
preferences.append(genPref([4,2,1,3], 0))
preferences.append(genPref([4,3,2,1], 0))
preferences.append(genPref([4,3,1,2], 0))

#similar to above, but order is a>c>d>b
#shows how ranking in non-schwartz set can't be computed
'''
preferences.append(genPref([1,2,3,4], 2))
preferences.append(genPref([1,2,4,3], 0))
preferences.append(genPref([1,3,2,4], 0))
preferences.append(genPref([1,3,4,2], 4))
preferences.append(genPref([1,4,2,3], 3))
preferences.append(genPref([1,4,3,2], 0))
preferences.append(genPref([2,1,4,3], 0))
preferences.append(genPref([2,1,3,4], 0))
preferences.append(genPref([2,3,4,1], 0))
preferences.append(genPref([2,3,1,4], 0))
preferences.append(genPref([2,4,3,1], 0))
preferences.append(genPref([2,4,1,3], 0))
preferences.append(genPref([3,1,2,4], 0))
preferences.append(genPref([3,1,4,2], 0))
preferences.append(genPref([3,2,1,4], 0))
preferences.append(genPref([3,2,4,1], 0))
preferences.append(genPref([3,4,1,2], 0))
preferences.append(genPref([3,4,2,1], 0))
preferences.append(genPref([4,1,3,2], 0))
preferences.append(genPref([4,1,2,3], 0))
preferences.append(genPref([4,2,3,1], 0))
preferences.append(genPref([4,2,1,3], 0))
preferences.append(genPref([4,3,2,1], 0))
preferences.append(genPref([4,3,1,2], 0))
'''

candMap = dict()
candMap[1] = 'alice'
candMap[2] = 'bob'
candMap[3] = 'charley'
candMap[4] = 'diane'

profile = Profile(candMap, preferences)
print(profile.getWmg())

mechanism1 = mechanism.MechanismSchulzeFW()
rankingFW = mechanism1.getRanking(profile)
print
print("Floyd-Warshall ranking:")
print(rankingFW)

mechanism2 = mechanism.MechanismSchulzeK()
rankingK = mechanism2.getRanking(profile)
print
print("Kosaraju ranking:")
print(rankingK)