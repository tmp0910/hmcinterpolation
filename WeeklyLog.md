# Weekly Log #

Work Log
  * Week 1
    * Began searching for a project paper
    * Found a potential paper called "Moving Gradients: A Path-Based Method for Plausible Image Interpolation" by Mahajan, Huang, Matusik, Ramamoorthi, and Belhumeur (sigraph 2009). A copy is archived at simonyanginfo.com/interpolation.pdf
    * Creation of the Project Proposal  (Schedule and Proposed Rubric)
    * Began in-depth reading of the main paper and supporting papers (graph cuts)
    * Individual write-ups of the overall algorithm
    * Some confusion remains for the energy minimization process (based on graph cuts; specifically, the graph structure)
    * Discussed what was written in the individual write-ups
    * Emailed the authors of the paper to see if we might be able to ask them questions or get sample testing images from them (didn't receive any responses)
  * Week 2
    * Worked on class presentation powerpoint
    * Did an in-class presentation
    * Cleared up confusion of path extension when transitioning to the next layer of the Gaussian pyramid
    * The laborious process of installing the PIL library for Python
    * Met with Prof. Z who recommended using a different energy minimization algorithm (a hill-climbing approach)
    * Decided to use C++ for increased computational speed
    * Initial Gaussian pyramid implementation written by Simon
    * Initial attempt to write the correspondence (local cost) portion of the energy function.
      * Uncertainties regarding the standard deviation computation
    * Further investigation of the local cost energy function
  * Week 3
    * Wrote the coherency energy component function
    * Wrote the code for the energy function
    * Coded for checking path validity
    * Completed initial path random assignment
    * Wrote the main interpolation loop
    * Code to run multiple iterations of path-finding to increase consistency of interpolated image result
    * Finished collecting/creating test images to test the program
    * Met with Prof. Z and clarified the definition of the standard deviation
    * Multiple Bug fixes, although this caused our algorithm to run slower
    * Final results presentation powerpoint
    * Wrote up the HTML writeup and consolidated our logs