"""
    Different result from MPS means that the fixed point of the equation

           |    |    |    |
        -- A -- A -- A -- A --         |    |    |    |
           |    |    |    |     ==  -- A -- A -- A -- A --
        -- O -- O -- O -- O --         |    |    |    |
           |    |    |    |

    can't rank-1 decompose exactly.

                        |
           |            C
        -- A --  ⇏     ⊗   (rank-1)
           |         -- A --
                        |

    It hints that there are topological difference between the case being open boundaries in two dimensions
    and the case being periodic boundary in one dimension and open boundary in the other dimension.
"""