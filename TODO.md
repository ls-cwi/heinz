High priority:

* Bug in multi-threading: cutCount = 0
* Branch on stars
* Break ties: smallest deg[i] - posDeg[i]
* If there is only one root node => use rooted separation...
* Encapsulate connected components separation in NodeCut class (argument x_values not needed)
* Unit test: solution is an entire connected component in the input graph
* If #separated cuts drops below threshold stop separating in the root node!

Low priority:

* More sophisticated branching rules (favor y-vars)
* Use rooted formulation in unrooted call backs (addLocal)
* Handle components inside callbacks (get rid of enumeration functionality)

Done:

* Only do min cut separation on nodes i that are not part of non-zero component containing the root
* Remove NodeCut::_root

Not done:

* Construct local h at every call! Slower...
