High priority:

* Get rid of distinction between run and construct
* Put bi/tri information in base class preprocessing rules
* Rooted formulation: identify S... and exploit
* Bug in multi-threading: cutCount = 0
* If there is only one root node => use rooted separation...
* Encapsulate connected components separation in NodeCut class (argument x_values not needed)
* Unit test: solution is an entire connected component in the input graph
* If #separated cuts drops below threshold stop separating in the root node!

New preprocessing rules:

* Look at biconnected components
* Diamond with two pos nodes at center and one corner with deg 2
* Positive deg 1 nodes, merge those (if they have a score below LB!)

Low priority:

* More sophisticated branching rules (favor y-vars)
* Use rooted formulation in unrooted call backs (addLocal)
* Handle components inside callbacks (get rid of enumeration functionality)

Done:

* Only do min cut separation on nodes i that are not part of non-zero component containing the root
* Remove NodeCut::_root

Not done:

* Construct local h at every call! Slower...
* Branch on stars... cplex does a better job
* Break ties: smallest deg[i] - posDeg[i]
