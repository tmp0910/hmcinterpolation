# Introduction #


Final Project Proposal
  * Sayuri Soejima
  * Simon Yang
  * Alex Yin


# Project Description #

Given two input images, the goal of the project is to implement a system that allows us to come up with a plausible interpolation of the two images. There are two components to our project that we must consider. The first portion is the path algorithm. The algorithm finds a plausible path for each pixel, and uses that for interpolation. In order to accomplish the interpolation process, we use energy minimization techniques to calculate transition points.

## General project work breakdown: ##
  * Path algorithm – 2 people
  * Energy minimization (Graph cuts) – 1 person



## Schedule (with weekly milestones): ##

11/21 (Saturday)
  * Get an understanding of the problem and the methods.
  * Put up our description write-ups (for both sections).
11/22 (Sunday)
  * Finish reading everyone else’s, send comments by email.
  * Also put up the changes on the wiki.
11/24 (Tuesday)
  * Have completed initial attempt at programming our parts.
  * Path implementation: path construction
  * Energy minimization: graph cuts
11/29 (Sunday)
  * Group meeting where we can have a work session/catch up
  * Path implementation: transition points
12/1 (Tuesday)
  * By this date, get a working implementation (without occlusion)
  * Energy minimization: Computing a local minimum
12/4 (Friday)
  * Group meeting for another work session/catch up
  * Path implementation: occlusion handling
12/7 (Monday)
  * Final group meeting for the HTML write-up.
  * Path implementation: Using gradients for interpolation
12/8 (Tuesday)
  * Due date for Final Project
  * Have everything ready to be turned in, including the HTML write-up.


# Bibliography: #


## Main paper ##

  * Moving gradients: a path-based method for plausible image interpolation
  * http://portal.acm.org/citation.cfm?id=1531348

## Additional papers ##

  * Fast Approximate Energy Minimization via Graph Cuts
  * http://www.cs.cornell.edu/rdz/papers/bvz-iccv99.pdf



  * Computing Visual Correspondence with Occlusions via Graph Cuts
  * http://www.cs.cornell.edu/~rdz/Papers/KZ-ICCV01-tr.pdf


# Suggested Grading Rubric: #
## 40 points: Path implementation ##
  * 5 points: Brief write-up from each member about the path implementation in the paper, to demonstrate our understanding. As part of this process, we will each read everyone’s write-up, in order to catch any misunderstandings before getting into heavy implementation.
  * 10 points: Path construction – determining possible paths from the input images.
  * 10 points: Computing transition points. This will involve energy-minimization, or min-cost.
  * 10 points: Occlusion handling.
  * 5 points: Use gradients instead of image intensities during interpolation, to improve the interpolation results (particularly with edges).

## 25 points: Energy function ##
  * 5 points: Brief write-up from each member about the energy function implemented in the paper, to demonstrate our understanding. As part of this process, we will each read everyone’s write-up, in order to catch any misunderstandings and to further our understanding.
  * 10 points: Graph cuts.
  * 10 points: Computing a local minimum.

## 10 points: HTML write-up describing the features we implemented ##
(75 total points)



# Details #

Add your content here.  Format your content with:
  * Text in **bold** or _italic_
  * Headings, paragraphs, and lists
  * Automatic links to other wiki pages