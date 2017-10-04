# Where's Croc Lab Assignment

+ moveInfo (A list of information for the move. This has two fiels. The first is a vector of numbers called 'moves', where you will enter the moves you want to make. You should enter two moves (so you can move to a neighboring waterhole and search). Valid moves are the numbers of a neighboring or current waterhole or '0' which means you will search your current waterhole for Croc. The second field is a list called 'mem' that you can use to store information you want to remember from turn to turn. ):
  - $moves
  - $mem

+ readings (A vector giving the salinity, phosphate and nitrogen reading from Croc sensors at his current location)
  - [1] salinity
  - [2] phosphate
  - [3] nitrogen

+ positions (A vector giving the positions of the two tourists and yourself. If a tourist has just been eaten by Croc that turn, the position will be multiplied by -1. If a tourist was eaten by Croc in a previous turn, then the position will be NA)
  - [0] tourist position
  - [1] tourist position
  - [2] our position

+ edges (a matrix giving the edges paths between waterholes (edges) present)

+ probs (a list of three matrices giving the mean and standard deviation of readings for salinity, phosphate and nitrogen respectively at each waterhole)
  - $salinity
  - $phosphate
  - $nitrogen
