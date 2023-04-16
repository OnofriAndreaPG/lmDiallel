## Changes in 'lmDiallel'

2021-06-22: First stable release, including both fixed and mixed effects diallel models.
2023-02-25: Corrected buglet in lm.diallel, which did not permit to use designs without blocking
2023-03-10: Corrected buglet in SCA(), which introduced an error for mating schemes without selfs
2023-03-24: Added some support for GRIFFING4 with missing crosses
2023-04-14: Solved some bugs that would prevent the display of correct results with ME diallel experiments with different missing crosses in different environments
2023-04-15: Solved a bug that caused some bugs resulting in different parameter estimates when the order of parentals was swapped in lm.diallel.
