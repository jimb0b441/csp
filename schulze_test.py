#to use for testing implementation of schulze

from profile import Profile
from preference import Preference
import mechanism 

#ranking 1: 2 @ a>b>c>d
wmgMap1 = dict()
wmgMap1[1] = dict()
wmgMap1[2] = dict()
wmgMap1[3] = dict()
wmgMap1[4] = dict()
wmgMap1[1][2] = 1
wmgMap1[1][3] = 1
wmgMap1[1][4] = 1
wmgMap1[2][1] = -1
wmgMap1[2][3] = 1
wmgMap1[2][4] = 1
wmgMap1[3][1] = -1
wmgMap1[3][2] = -1
wmgMap1[3][4] = 1
wmgMap1[4][1] = -1
wmgMap1[4][2] = -1
wmgMap1[4][3] = -1
preference1 = Preference(wmgMap1,2)

#ranking 2: 3 @ b>d>c>a
wmgMap2 = dict()
wmgMap2[1] = dict()
wmgMap2[2] = dict()
wmgMap2[3] = dict()
wmgMap2[4] = dict()
wmgMap2[1][2] = -1
wmgMap2[1][3] = -1
wmgMap2[1][4] = -1
wmgMap2[2][1] = 1
wmgMap2[2][3] = 1
wmgMap2[2][4] = 1
wmgMap2[3][1] = 1
wmgMap2[3][2] = -1
wmgMap2[3][4] = -1
wmgMap2[4][1] = 1
wmgMap2[4][2] = -1
wmgMap2[4][3] = 1
preference2 = Preference(wmgMap2,3)

#ranking 3: 5 @ C>D>B>A
wmgMap3 = dict()
wmgMap3[1] = dict()
wmgMap3[2] = dict()
wmgMap3[3] = dict()
wmgMap3[4] = dict()
wmgMap3[1][2] = -1
wmgMap3[1][3] = -1
wmgMap3[1][4] = -1
wmgMap3[2][1] = 1
wmgMap3[2][3] = -1
wmgMap3[2][4] = -1
wmgMap3[3][1] = 1
wmgMap3[3][2] = 1
wmgMap3[3][4] = 1
wmgMap3[4][1] = 1
wmgMap3[4][2] = 1
wmgMap3[4][3] = -1
preference3 = Preference(wmgMap3,5)

#ranking 4: 6 @ a>c>b>d
wmgMap4 = dict()
wmgMap4[1] = dict()
wmgMap4[2] = dict()
wmgMap4[3] = dict()
wmgMap4[4] = dict()
wmgMap4[1][2] = 1
wmgMap4[1][3] = 1
wmgMap4[1][4] = 1
wmgMap4[2][1] = -1
wmgMap4[2][3] = -1
wmgMap4[2][4] = 1
wmgMap4[3][1] = -1
wmgMap4[3][2] = 1
wmgMap4[3][4] = 1
wmgMap4[4][1] = -1
wmgMap4[4][2] = -1
wmgMap4[4][3] = -1
preference4 = Preference(wmgMap4,6)

preferences =[preference1, preference2, preference3, preference4]
candMap = dict()
candMap[1] = 'alice'
candMap[2] = 'bob'
candMap[3] = 'charley'
candMap[4] = 'diane'

profile = Profile(candMap, preferences)
print(profile.getWmg())
mechanism = mechanism.MechanismSchulzeK()
print(mechanism.getRanking(profile))
