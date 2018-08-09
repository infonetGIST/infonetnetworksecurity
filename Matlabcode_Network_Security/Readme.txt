This pacakge includes Compromise Detection and Attack Compensation algorithms in relay-assisted wireless multiple access networks with two-phase transmission protocols.

The following codes can produce some of simulation results in our paper titled as Disabling Pollution Attacks on Cooperative Wireless Multiple Access Networks.

******************** main codes ********************
0. H_Generator               :  It generates a structure of network coding
1. Contradiction_analysis   :   It generates contradiction probabilities of compromised relays and usual relays
2. Observation_window     :  It generates observation window (N) to satisfy Neyman-pearson critieron, i.e., false alarm probability
3. Detection_evaluation    :   It evaluates perfomance of the compromise detection with determined observation window 
4. Performance_evaluation :   It evalutes performance of proposed compromise detection and attack compensation methods

********************     Data     ********************
H_matrix : Generating matrix and parity-check matrix obtained in '0. H_Generator' are saved.

********************   functions  ********************
NC_decoder :  functions of network decoding 
MP_analysis  :  functions of message-passing decoding
