# quadrilateral-element-4to9nodes-fortran

== Analyzing Quadrilateral elements with 4-9 nodes using FORTRAN90 ==

==== Developed by Mohammad Sajednia ====

Inputs:

1. nodes.txt
Column 1: node number
Column 2: X coordinate of the node
Column 3: Y coordinate of the node

2. elements.txt
Column 1: the number of the element
Column 2-10: the number of nodes

Note 1: If the element doesn't include 9 nodes, fill the remaining columns with zeros.

Note 2: the order of numbering nodes is like following:

3===6===2

=---------=

7----9----5

=---------=

4===8===1

3. boundaries.txt
Column 1: the number of the constrained node
Column 2: the direction of the constraint (1 for X, 2 for Y)
Column 3: the value of the constraint (0 if the constraint is a support and any value for support settlements)

4. nloads.txt
Column 1: the number of the node
Column 2: the magnitude of the nodal force
Column 3: the angle between the force vector and the positive direction of the X axis

5. tractions.txt
Column 1: the number of nodes on the edge (2or3)
Column 2: the number of the starting node
Column 3: the number of the middle node (0 if the edge includes only 2 nodes)
Column 4: the number of the ending node
Column 5: the magnitude of the traction on the starting node
Column 6: the magnitude of the traction on the middle node (0 if the edge includes only 2 nodes)
Column 7: the magnitude of the traction on the ending node
Column 8: the angle between the traction vector and the positive direction of the X axis

6. vloads.txt
Column 1: the number of the element that includes the body force
Column 2: the magnitude of the body force
Column 3: the angle between the body force vector and the positive direction of the X axis
