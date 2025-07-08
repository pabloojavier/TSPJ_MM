# Branch-and-cut algorithms for the traveling salesman problem with job-times
date: "July, 2025"

## Abstract:

The traveling salesman problem with job times (TSPJ) involves a traveler visiting a set of vertices, ensuring one visit to each of them while he initiates a job. The time of each job depends on the vertex where it is performed. Once started, the traveler moves to the next vertex, and the job continues autonomously. The aim is to minimize the maximum completion time. Since the problem is NP-hard, existing mixed-integer linear programming (MILP) models are unable to solve large instances efficiently. Therefore, this paper proposes to enhance the existing MILP model and introduces a new MILP model for the TSPJ, which incorporates valid lower and upper bounds to strengthen them. Moreover, we propose two branch-and-cut (B\&C) algorithms based on the improved existing model and the new one. These algorithms integrate strengthened exponential-size formulations that explicitly incorporate subtour elimination constraints and blossom inequalities. B\&C algorithms are tested on instances from the literature, comprising four sets of instances with sizes ranging from 17 to 1200 vertices. Computational results show that the proposed B\&C algorithms outperform the state-of-the-art MILP models in all instances. Since no formulation achieves optimality within the given time limit for large instances, only three medium instances with up to 454 vertices have reached optimality. Consequently, we propose a fifth set of instances, ranging from 100 to 390 vertices, to further assess performance limits. B\&C algorithms demonstrate improved performance with lower gap values in all instances, and faster computing times while optimally solving instances with up to 386 vertices.


## How to run
* All codes are in 'Codigos' folder.
* Just run 'MM_po.py' file to run just one instance
* To run all instance from one size, run 'multi.py'