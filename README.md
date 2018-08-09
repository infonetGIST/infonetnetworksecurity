# infonet network security

This pacakge includes Compromise Detection and Attack Compensation algorithms in relay-assisted wireless multiple access networks with two-phase transmission protocols.

The following codes can produce some of simulation results in our paper titled as ***Disabling Pollution Attacks on Cooperative Wireless Multiple Access Networks***.

# main codes 

0. H_Generator: generates a generator matrix and a parity check matrix given the structure of cooperative wireless multiple access network.
1. Contradiction_analysis: generates contradiction probabilities of compromised relays and usual relays.
2. Observation_window: generates observation window (*N*) to satisfy Neyman-Pearson criterion, i.e., false alarm probability.
3. Dection_evaluation: evalutes performance of proposed comporomise detection and attack compensation methods.

# data 

 H_matrix: Generator matrix and parity-check matrix obtained in 'H_Generator' are stored.

# functions

NC_decoder: functions of network deocding.

MP_analysis: functions of message-passing decoding for analysis.

# Copyright

The original package is available at http://infonet.gist.ac.kr/

COPYRIGHT (c) 2018 Heung-No Lee, and Woong-Bi Lee.

E-mail: <heungno@gist.ac.kr>, woongbi.lee@gmail.com

This package is distributed under the terms of the GNU General Public License 3.0. http://www.gnu.org/copyleft/gpl.html

Permission to use, copy, modify, and distribute this software for any purpose without fee is hereby granted, provided that this entire notice is included in all copies of any software which is or includes a copy or modification of this software and in all copies of the supporting documentation for such software. This software is being provided "as is", without any express or implied warranty. In particular, the authors do not make any representation or warranty of any kind concerning the merchantability of this software or its fitness for any particular purpose.