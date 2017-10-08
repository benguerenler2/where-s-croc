# Authors:
# Arianna Delsante - 940929-T300
# Bengü Erenler - 940519-T520
# Ankur Shukla - 950803-4734
# Diego Castillo - 911206-T438

###############################################
#############Inference Related Code############
###############################################

# Return the previous state for the probability of the Croc being at
# any of the available waterholes
getPrevState=function(prevState, numOfWaterHoles, positions) {
  if(isTRUE(is.null(prevState))) {
    backpackerA = positions[1]
    backpackerB = positions[2]
    if (backpackerA == backpackerB) {
      matrix = matrix(1/(numOfWaterHoles - 1), nrow=1, ncol=numOfWaterHoles)
      matrix[backpackerA] = 0
      return (matrix)
    } else {
      matrix = matrix(1/(numOfWaterHoles - 2), nrow=1, ncol=numOfWaterHoles)
      matrix[backpackerA] = 0
      matrix[backpackerB] = 0
      return (matrix)
    }
  }
  return (prevState)
}

# Create the transition matrix for each available waterhole
createTransitions = function(edges, numOfWaterHoles) {
  matrix = matrix(0, nrow=numOfWaterHoles, ncol=numOfWaterHoles)
  for (waterhole in 1:numOfWaterHoles){
    options = getOptions(waterhole, edges)
    prob = 1/length(options)
    for (option in options) {
      matrix[option, waterhole] = prob
    }
  }
  return (matrix)
}

# Create an observation matrix given readings of the current position
# of the Croc and a normal distribution of the readings throughout all
# of the waterholes
createObservations=function(readings, probs, numOfWaterHoles) {
  matrix = matrix(0, nrow=numOfWaterHoles, ncol=numOfWaterHoles)
  for (waterhole in 1:numOfWaterHoles) {
    salinityDNorm = dnorm(readings[1], probs$salinity[waterhole,1], probs$salinity[waterhole,2])
    phosphateDNorm = dnorm(readings[2], probs$phosphate[waterhole,1], probs$phosphate[waterhole,2])
    nitrogenDNorm = dnorm(readings[3], probs$nitrogen[waterhole,1], probs$nitrogen[waterhole,2])
    matrix[waterhole, waterhole] = salinityDNorm * phosphateDNorm * nitrogenDNorm
  }
  return (matrix)
}

# Returns a normalized state
normalizeState=function(state, numOfWaterHoles) {
  stateTotal = sum(state)
  normalizedState = matrix(0, nrow=1, ncol=numOfWaterHoles)
  for (i in 1:numOfWaterHoles) {
    normalizedState[i] = state[i]/stateTotal
  }
  return (normalizedState)
}

# Use the forward algorithm to provide a distribution over the system
getInferedState=function(prevState, readings, edges, probs, numOfWaterHoles, positions) {
  # Initialize previous state, transition, and observation matrices
  prevState = getPrevState(prevState, numOfWaterHoles, positions)
  transitions = createTransitions(edges, numOfWaterHoles)
  observations = createObservations(readings, probs, numOfWaterHoles)

  # Compute next state
  return (prevState %*% (transitions %*% observations))
}

# Return board state given that we know a backpacker was eaten by the Croc
# on this turn
getEatenBackpackerState=function(positions, numOfWaterHoles) {
  goal = abs(positions[positions[] < 0][1])
  state = matrix(0, nrow=1, ncol=numOfWaterHoles)
  state[goal] = 1
  return (state)
}

# Return true if a backpacker has been eaten in this turn
hasBackpackerBeenEaten=function(positions) {
  allNumbers = isTRUE(all(positions == floor(positions)))
  hasBeenEaten = length(positions[positions[] < 0]) > 0
  return (allNumbers && hasBeenEaten)
}

# Given backpacker's position, return true if backpacker is alive, false otherwise
isBackpackerAlive=function(position) {
  return (isTRUE(position == floor(position)) && position > 0)
}

# Given a path, return the best next move towards goal
generateNextMove=function(path) {
  if(isTRUE(length(path) == 1)) {
    return (path[[1]])
  } else {
    return (path[[2]])
  }
}

hmmWC=function(moveInfo, readings, positions, edges, probs) {
  numOfWaterHoles = dim(probs$phosphate)[1]
  prevState = moveInfo$mem$prevState
  currState = NULL

  # Use different strategies to determine the current state probability distribution
  if(hasBackpackerBeenEaten(positions)){
    currState = getEatenBackpackerState(positions, numOfWaterHoles)
  } else {
    currState = getInferedState(prevState, readings, edges, probs, numOfWaterHoles, positions)
  }

  # Search
  currState = normalizeState(currState, numOfWaterHoles)
  from = positions[3]
  goal = which.max(currState)
  path = bestFirstSearch(from, goal, edges, getPoints(), currState)

  # Generate next move, pass next turn previous state
  moveInfo$mem$prevState = currState
  moveInfo$moves = c(generateNextMove(path), 0)
  return (moveInfo)
}

averageTest <- function(tests, showCroc = FALSE, pause = 0, doPlot = FALSE){
  sum = 0
  for (i in 1:tests) {
    sum=sum+runWheresCroc(hmmWC, showCroc, pause, doPlot)
    if(i%%10==0){
      print(i)
      print(sum/i)
    }
  }
  print(sum/i)
  return(0)
}

###############################################
##############Search Related Code##############
###############################################

# A priority queue which allows to insert elements
# and order them by priority
# Source: http://rosettacode.org/wiki/Priority_queue#R
PriorityQueue <- function() {
  queueKeys <<- queueValues <<- NULL
  insert <- function(key, value) {
    # If node already exists on queue, and this new addition is better,
    # delete previous one and insert this new one instead
    index = getValueIndex(value)
    if(length(index) > 0) {
      if(isTRUE(key < queueKeys[[index]])) {
        queueKeys <<- queueKeys[-index]
        queueValues <<- queueValues[-index]
      } else {
        # Ignore it, we already have a cheaper path
        return (-1)
      }
    }

    # Insert new value in queue
    temp <- c(queueKeys, key)
    ord <- order(temp)
    queueKeys <<- temp[ord]
    queueValues <<- c(queueValues, list(value))[ord]
  }
  pop <- function() {
    head <- queueValues[[1]]
    queueValues <<- queueValues[-1]
    queueKeys <<- queueKeys[-1]
    return (head)
  }
  empty <- function() length(queueKeys) == 0
  getValueIndex <- function(value) which(queueValues %in% list(value) == TRUE)
  list(insert = insert, pop = pop, empty = empty)
}

# A simple lists which allows to insert elements on it
# and verity if a particular element exists or not
List <- function() {
  listValues <<- NULL
  insert <- function(value) listValues <<- c(listValues, list(value))
  exists <- function(value) isTRUE(which(listValues %in% list(value) == TRUE) > 0)
  list(insert = insert, exists = exists)
}

# Return the neighbors of a particular node
getNeighbors=function(node, edges) {
  neighbors = getOptions(node, edges)
  # A node is not a neighbor of itself
  return (neighbors[neighbors[] != node])
}

# Return true if node is goal, false otherwise
isGoal=function(node, goal) {
  return (node == goal)
}

# Transform a vector representation of a node to a string
transformNodeToString=function(node) {
  return (paste(node))
}

# Transform a string representation of a node to a vector
transformStringToNode=function(nodeAsString) {
  return (c(as.integer(nodeAsString)))
}

# Returns the path from an initial position to the goal position given
# the path visited by the algorithm
generatePath=function(from, to, path) {
  goal = transformNodeToString(from)
  curr = transformNodeToString(to)

  # Build path visited by traversing the path variable
  # from goal to initial position (in reverse order)
  vectors = list(c(to))
  while (curr != goal) {
    node = transformStringToNode(path[[curr]])
    vectors = c(vectors, list(node))
    curr = path[[curr]]
  }

  # Return path from initial position to goal
  return (rev(vectors))
}

# Add a node to the path
addNodeToPath=function(path, from, to) {
  path[transformNodeToString(to)] = transformNodeToString(from)
  return (path)
}

# Return the heuristic value of going from a node to another one
getHeuristicValue=function(from, node, neighbor, path, state) {
  # Compute edge cost
  neighborPath = addNodeToPath(path, node, neighbor)
  edgeCost = length(generatePath(from, neighbor, neighborPath)) - 1

  # Favor nodes which have a higher likelihood for the Croc being there
  return (edgeCost - state[neighbor])
}

# Find goal though a best-first search
bestFirstSearch=function(from, goal, edges, locations, state) {
  # Initialize visited, frontier, and path lists
  visited = List()
  frontier = PriorityQueue()
  path = list()

  # Put the starting location on the frontier (heuristic 0 is fine)
  frontier$insert(0, from)

  while (!frontier$empty()) {
    # Get node with the least heuristic on the frontier
    node = frontier$pop()

    # Return the visited path + current node as path to goal
    if(isGoal(node, goal)) {
      return (generatePath(from, goal, path))
    }

    neighbors = getNeighbors(node, edges)
    for (neighbor in neighbors) {
      # Only search neighbors which haven't already being visited
      if(visited$exists(neighbor)) {
        next
      } else {
        heuristic = getHeuristicValue(from, node, neighbor, path, state)
        inserted = frontier$insert(heuristic, neighbor)

        # Add neighbor to path only if it was inserted in the frontier
        wasInserted = length(inserted) != 1 || inserted[[1]][1] != -1
        if (isTRUE(wasInserted)) {
          path = addNodeToPath(path, node, neighbor)
        }
      }
    }

    # Keep track of visited nodes
    visited$insert(node)
  }
}

#' @export
randomWC=function(moveInfo,readings,positions,edges,probs) {
  moveInfo$moves=c(sample(getOptions(positions[3],edges),1),0)
  return(moveInfo)
}

#' @export
manualWC=function(moveInfo,readings,positions,edges,probs) {
  options=getOptions(positions[3],edges)
  print("Move 1 options (plus 0 for search):")
  print(options)
  mv1=readline("Move 1: ")
  if (mv1=="q") {stop()}
  if (!mv1 %in% options && mv1 != 0) {
    warning ("Invalid move. Search ('0') specified.")
    mv1=0
  }
  if (mv1!=0) {
    options=getOptions(mv1,edges)
  }
  print("Move 2 options (plus 0 for search):")
  print(options)
  mv2=readline("Move 2: ")
  if (mv2=="q") {stop()}
  if (!mv1 %in% options && mv1 != 0) {
    warning ("Invalid move. Search ('0') specified.")
    mv2=0
  }
  moveInfo$moves=c(mv1,mv2)
  return(moveInfo)
}

#' Run Where's Croc
#'
#' Runs the Where's Croc game. In this game, you are a ranger in an Australian national park.
#' This park consists of a number of waterholes, some of which are connected to each other.
#' There is a crocodile in the park called 'Croc'. Croc has been fitted with sensors that record
#' the salinity, phosphate and nitrogen levels in the water where he is swimming. He was also
#' fitted with a sensor that records his position, but that has broken.
#' Your task is to find Croc using the available information. To aid in this you have information
#' about the probability distributions for different salinity, phosphate and nitrogen levels in
#' different waterholes.
#' There are also two tourists in the park. Both the tourists and Croc walk randomly, each turn
#' moving to one of the neighboring waterholes from where they are or staying still. All moves
#' are equally likely.
#' If Croc and a tourist end up on the same waterhole, Croc will eat the tourist. If you search
#' the waterhole you are on when Croc is there, you have found Croc and win the game.
#' Your score is the number of turns it takes to find Croc.
#' To play manually pass manualWC
#' as the makeMoves function and enter the appropriate numbers to make moves.
#' @param makeMoves Your function that takes five arguments: (1) A list of information for the move.
#' This has two fiels. The first is a vector of numbers called 'moves', where you will enter
#' the moves you want to make. You should
#' enter two moves (so you can move to a neighboring waterhole and search). Valid moves are the
#' numbers of a neighboring or current waterhole or '0' which means you will search your current
#' waterhole for Croc. The second field is a list called
#' 'mem' that you can use to store information you want to remember from turn to turn. (2) A
#' vector giving the salinity, phosphate and nitrogen reading from Croc sensors at his current
#' location. (3) A vector giving the positions of the two tourists and yourself. If a tourist
#' has just been eaten by Croc that turn, the position will be multiplied by -1. If a tourist
#' was eaten by Croc in a previous turn, then the position will be NA. (4) a matrix giving the
#' edges paths between waterholes (edges) present. (5) a list of three matrices giving the mean
#' and standard deviation of readings for salinity, phosphate and nitrogen respectively
#' at each waterhole.
#' Your function should return the first argument passed with an updated moves vector
#' and any changes to the 'mem' field you wish to access later on.
#' @param showCroc A Boolean value specifying whether you want Croc to be shown on the gameboard.
#' Note that you are not permitted to use this visual information when you are scored.
#' @param pause The pause period between moves. Ignore this.
#' @return A string describing the outcome of the game.
#' @export
runWheresCroc=function(makeMoves,showCroc=F,pause=1, doPlot=TRUE) {
  positions=sample(1:40,4) # Croc, BP1, BP2, Player
  points=getPoints()
  edges=getEdges()
  probs=getProbs()
  move=0
  moveInfo=list(moves=c(),mem=list())
  while (!is.na(positions[1])) {
    move=move+1
    positions[1]=sample(getOptions(positions[1],edges),1)
    if (!is.na(positions[2])&&positions[2]>0) {
      positions[2]=sample(getOptions(positions[2],edges),1)
    } else if (!is.na(positions[2]) && positions[2]<0) {
      positions[2]=NA
    }
    if (!is.na(positions[3])&&positions[3]>0) {
      positions[3]=sample(getOptions(positions[3],edges),1)
    } else if (!is.na(positions[3]) && positions[3]<0) {
      positions[3]=NA
    }
    if (!is.na(positions[2]) && positions[2]==positions[1]) {
      positions[2]=-positions[2]
    }
    if (!is.na(positions[3]) && positions[3]==positions[1]) {
      positions[3]=-positions[3]
    }

    if(doPlot) {
      plotGameboard(points,edges,move,positions,showCroc)
    }

    Sys.sleep(pause)

    readings=getReadings(positions[1],probs)
    moveInfo=makeMoves(moveInfo,readings,positions[2:4],edges,probs)
    if (length(moveInfo$moves)!=2) {
      stop("Error! Passed makeMoves function should return a vector of two elements.")
    }
    for (m in moveInfo$moves) {
      if (m==0) {
        if (positions[1]==positions[4]) {
          print(paste("Congratulations! You got croc at move ",move,".",sep=""))
          return (move)
        }
      } else {
        if (m%in%getOptions(positions[4],edges)) {
          positions[4]=m
        } else {
          warning("Invalid move.")
        }
      }
    }
  }
}
#' @export
getPoints=function() {
  points=matrix(c(1,1),ncol=2)
  points=rbind(points,c(1,7))
  points=rbind(points,c(1,17))
  points=rbind(points,c(2,3))
  points=rbind(points,c(2,12))
  points=rbind(points,c(3,2))
  points=rbind(points,c(3,19))
  points=rbind(points,c(4,7))
  points=rbind(points,c(4,11))
  points=rbind(points,c(5,5))
  points=rbind(points,c(5,15))
  points=rbind(points,c(6,1))
  points=rbind(points,c(6,20))
  points=rbind(points,c(7,6))
  points=rbind(points,c(7,11))
  points=rbind(points,c(8,2))
  points=rbind(points,c(8,14))
  points=rbind(points,c(8,18))
  points=rbind(points,c(9,6))
  points=rbind(points,c(10,10))
  points=rbind(points,c(10,18))
  points=rbind(points,c(11,1))
  points=rbind(points,c(11,12))
  points=rbind(points,c(12,6))
  points=rbind(points,c(12,12))
  points=rbind(points,c(13,16))
  points=rbind(points,c(14,4))
  points=rbind(points,c(14,12))
  points=rbind(points,c(14,20))
  points=rbind(points,c(15,3))
  points=rbind(points,c(15,8))
  points=rbind(points,c(15,17))
  points=rbind(points,c(16,14))
  points=rbind(points,c(17,3))
  points=rbind(points,c(17,18))
  points=rbind(points,c(18,10))
  points=rbind(points,c(19,13))
  points=rbind(points,c(20,2))
  points=rbind(points,c(20,6))
  points=rbind(points,c(20,19))
  return (points)
}

#' @export
getEdges=function() {
  edges=matrix(c(1,2),ncol=2)
  edges=rbind(edges,c(1,4))
  edges=rbind(edges,c(1,6))
  edges=rbind(edges,c(2,4))
  edges=rbind(edges,c(2,5))
  edges=rbind(edges,c(3,5))
  edges=rbind(edges,c(3,7))
  edges=rbind(edges,c(4,6))
  edges=rbind(edges,c(4,8))
  edges=rbind(edges,c(5,7))
  edges=rbind(edges,c(5,9))
  edges=rbind(edges,c(6,12))
  edges=rbind(edges,c(7,11))
  edges=rbind(edges,c(7,13))
  edges=rbind(edges,c(8,9))
  edges=rbind(edges,c(8,10))
  edges=rbind(edges,c(9,11))
  edges=rbind(edges,c(10,12))
  edges=rbind(edges,c(10,14))
  edges=rbind(edges,c(11,13))
  edges=rbind(edges,c(11,15))
  edges=rbind(edges,c(12,16))
  edges=rbind(edges,c(13,18))
  edges=rbind(edges,c(14,15))
  edges=rbind(edges,c(14,16))
  edges=rbind(edges,c(15,17))
  edges=rbind(edges,c(16,19))
  edges=rbind(edges,c(16,22))
  edges=rbind(edges,c(17,18))
  edges=rbind(edges,c(17,19))
  edges=rbind(edges,c(17,20))
  edges=rbind(edges,c(18,21))
  edges=rbind(edges,c(19,20))
  edges=rbind(edges,c(19,22))
  edges=rbind(edges,c(20,23))
  edges=rbind(edges,c(21,23))
  edges=rbind(edges,c(21,29))
  edges=rbind(edges,c(22,24))
  edges=rbind(edges,c(22,27))
  edges=rbind(edges,c(23,24))
  edges=rbind(edges,c(23,25))
  edges=rbind(edges,c(24,25))
  edges=rbind(edges,c(24,27))
  edges=rbind(edges,c(25,26))
  edges=rbind(edges,c(25,27))
  edges=rbind(edges,c(25,28))
  edges=rbind(edges,c(26,28))
  edges=rbind(edges,c(26,29))
  edges=rbind(edges,c(27,30))
  edges=rbind(edges,c(27,31))
  edges=rbind(edges,c(28,31))
  edges=rbind(edges,c(28,32))
  edges=rbind(edges,c(29,32))
  edges=rbind(edges,c(29,35))
  edges=rbind(edges,c(30,31))
  edges=rbind(edges,c(30,34))
  edges=rbind(edges,c(31,33))
  edges=rbind(edges,c(31,34))
  edges=rbind(edges,c(32,33))
  edges=rbind(edges,c(32,35))
  edges=rbind(edges,c(33,35))
  edges=rbind(edges,c(33,36))
  edges=rbind(edges,c(33,37))
  edges=rbind(edges,c(34,36))
  edges=rbind(edges,c(34,38))
  edges=rbind(edges,c(35,40))
  edges=rbind(edges,c(36,37))
  edges=rbind(edges,c(36,39))
  edges=rbind(edges,c(37,39))
  edges=rbind(edges,c(37,40))
  edges=rbind(edges,c(38,39))

  return (edges)
}

#' @export
getProbs=function(){
  salinity=cbind(runif(40,100,200),runif(40,5,30))
  phosphate=cbind(runif(40,100,200),runif(40,5,30))
  nitrogen=cbind(runif(40,100,200),runif(40,5,30))
  list(salinity=salinity,phosphate=phosphate,nitrogen=nitrogen)
}

#' @export
getReadings=function(point,probs){
  c(
    rnorm(1,probs$salinity[as.numeric(point),1],probs$salinity[as.numeric(point),2]),
    rnorm(1,probs$phosphate[as.numeric(point),1],probs$phosphate[as.numeric(point),2]),
    rnorm(1,probs$nitrogen[as.numeric(point),1],probs$nitrogen[as.numeric(point),2])
  )
}


#' @export
plotGameboard=function(points,edges,move,positions,showCroc) {
  plot(points,pch=18,col="blue",cex=2,xlab="X",ylab="Y",main=paste("Where's Croc - Move",move))
  xFrom=points[edges[,1],1]
  yFrom=points[edges[,1],2]
  xTo=points[edges[,2],1]
  yTo=points[edges[,2],2]
  segments(xFrom,yFrom,xTo,yTo)
  for (bp in 2:3)
    if (!is.na(positions[bp])) {
      if (positions[bp]>0) {
        points(points[as.numeric(positions[bp]),1],points[as.numeric(positions[bp]),2],col="orange",pch=17,cex=4)
      } else {
        points(points[-as.numeric(positions[bp]),1],points[-as.numeric(positions[bp]),2],col="red",pch=17,cex=4)
      }
    }
  points(points[as.numeric(positions[4]),1],points[as.numeric(positions[4]),2],col="green",pch=15,cex=4)
  if (showCroc) {
    points(points[as.numeric(positions[1]),1],points[as.numeric(positions[1]),2],col="red",pch=15,cex=4)
  }
  text(points[,1]+.4, points[,2], labels=as.character(1:40))
}

#' @export
getOptions=function(point,edges) {
  c(edges[which(edges[,1]==point),2],edges[which(edges[,2]==point),1],point)
}
