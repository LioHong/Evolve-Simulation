# Evolve: An Artificial Life Simulator

## Contents

### Maps

1. big_bang

2. big_bang1

3. bot_battle1

4. bot_battle2

5. bot_battle3

6. bot_battle4

7. bot_battle5

8. default

9. default1

10. draw

11. eight_strains

12. eight_strains1

13. gonad

14. hanoi

15. intestines

16. intestines1

17. just_ducks

18. maze

19. snowflake

### Scripts

1. attacker

2. blueprint2

3. cmove

4. contest

5. diva

6. draw

7. duck

8. duck_killer

9. duck_killer2

10. gonad

11. greedy_seed

12. hanoi

13. hanoi15

14. iamthestar

15. iamthestar2

16. knut

17. knut2

18. knut3

19. seed

20. sex_seed

21. shell_diva

22. sperm

23. sqr

24. sqr2

25. sweeper

26. weasel

27. worm

28. xy_seed

### Example Scripts ('kforth_example_')

1. add

2. bigints

3. cube

4. factorial

5. factorial_loop

6. foobar

7. goo

8. infloop

9. junk

10. junk2

11. mul

12. mycall

13. myfirstprog

14. myif

15. pack

16. repeat

17. saveregs

18. simple

19. small

20. table

## Documentation of Scripts

### attacker

Moves toward stuff and eats it.

main:

reproduce call

Find nearest organic matter. Then try to eat it and towards it.

Loop

reproduce:

compute 1/4 of our energy, and store in R0

Make 1st spore to the square on our left.

Put 2nd spore at the same spot (fertilizing it)

Move in opposite direction from nearest spore

### blueprint2

main:

    NewCell call

    LiveLife call

LiveLife:

    EatForEver

    BeBone

    BeBrain

    BeFoot

All Be---- call Expand.

Expand:

    Blue_read

Blue_read:

    Blueprint

Blueprint:

    Diagram of organism shape

### cmove

This creature demonstrates CMOVE. It is based on the documentation for the CMOVE instruction.
This creature will morph from the shape shown in STEP 1, to the creature shown in STEP 4.

To syncronize activities between cells, the MOOD register in Cell 1 will be used.

STEP 1: Firstly, build the creature to look like this:

    +---+---+
    | 1 | 2 |
    +---+---+---+
        | 3 | 4 |
        +---+---+
        | 5 |
        +---+

STEP 2: cell 1 executes:

    0 1 CMOVE

Giving,

        +---+
        | 2 |
    +---+---+---+
    | 1 | 3 | 4 |
    +---+---+---+
        | 5 |
        +---+

STEP 3: cell 4 waits for the MOOD of cell 1 to become the value '4', then executes:

    0 1 CMOVE

Giving,

        +---+
        | 2 |
    +---+---+
    | 1 | 3 |
    +---+---+---+
        | 5 | 4 |
        +---+---+

STEP 4: cell 5 waits for the MOOD of cell 1 to become the value '5', then executes:

    1 -1 CMOVE

Giving,

        +---+
        | 2 |
    +---+---+---+
    | 1 | 3 | 5 |
    +---+---+---+
            | 4 |
            +---+

Finally, the the creature spins forever (using ROTATE).

Subfunctions:

- Cell1

- Cell2

- Cell3

- Cell4

- Cell5

- forever

- wait

- short_wait

### contest

This is the Contest Creature

CONTEST: Human Designed KFORTH versus Evolve KFORTH

Taken from the day-0068.evolve simulation file

Organism ID: 14,707,290,689
Location: (402,335)
Step: 222,544,721

CONTEST RULES:

1. Place your organism in a universe with this one.
2. Your ogranism is strain 1 with 10,000 units of energy
3. This organism is strain 2 with 10,000 units of energy
4. The universe is 700 x 600 with just the oval barrier
5. Mutations use the default settings

GOAL:
To cause strain 2 to go extinct. No time limit is
given.

### diva

Only goal is to move itself to the center

main:

    X-movement

    Y-movement

    Length to go

    Length to stop (4 distance)

    Distance

    Duck seen

    doturns call

    forever call

lookforduck:

    ; Look for nearest cell and store in (R0, R1) and call approach

doturns

    ; walk R4 moves

    lookforduck call

    ; repeat

walkline

    ; (length of move -- ) 

    ; go one cell

    ; look for duck

    ; You never know :-) - Try and eat

    ; decrease input one

    ; <>0 repeat

    ; pop input

forever

approach

    ; go one cell

    ; decrease input one

    ; <>0 repeat

    ; pop line

### draw

Do some drawing

main:

    ; X-movement

    ; Y-movement

    ; Tricks to go

    ; Length to go 

    doturns

    forever

forever

line

turnright

turnleft

move

draw

### duck

; Simple look/move/eat creature

main:

    ; Eat while something edible is 2 units away

    ; Find the farthest thing, and store direction in (R0, R1)

    ; Move in this direction, until stopped

    ; Find the farthest barrier (within a 10 units radius) and store direction in (R0, R1)

    ; Move in this direction, until stopped

### duck_killer

This creature is a 9 cell organism (WITHOUT healing ability). The leaf cells just eat like crazy, while the brain cell LOOK's and MOVE's toward the closest living thing.

      +-------+-------+-------+
      |       |       |       |
      | leaf5 | leaf1 | leaf7 |
      |       |       |       |
      +-------+-------+-------+
      |       |       |       |
      | leaf3 | brain | leaf4 |
      |       |       |       |
      +-------+-------+-------+
      |       |       |       |
      | leaf8 | leaf2 | leaf6 |
      |       |       |       |
      +-------+-------+-------+

The brain looks in all 8 directions (NW, N, NE, E, SE, S, SW, W) and moves in the direction that reports the smallest distance to a 'CELL'. If none of the 8 directions report a "hit" then use the "hunt" routine to move in a predefined pattern.

main:

    ; Try to grow and `leaf` in all directions

    ; initialize 1st direction for hunt mode.

    ; (R7, R8) is the most recent successful look vector

brain:

    ; Look in all 8 directions and move toward the one that has the closest thing. If none of them report a "hit" then call the 'hunt' routine to move in a predefined pattern.

direction_array:
    ; Index this array to translate a number (1 thru 8) into an (x,y) vector that can be used with OMOVE.

leaf:

    ; EAT routine. Try to maximize the number of EAT instructions versus the non-eating instructions.

hunt:

    When all the 8 LOOK operations return no hits, then we must move along some kind of hunting pattern.

    First, If (R7, R8) is non-zero we move in that direction, otherwise we move in the direction specified by R5. If the move using R5 fails we incrment R5 to the next direction.

    Registers used:
        R5 - contains the current direction index
        (R7, R8) most recently sucessful vector

### duck_killer2

This creature is a 12 cell organism without healing ability. The leaf cells just eat like crazy, while the brain cell LOOK's and MOVE's toward the closest living thing.

      +-------+-------+-------+-------+
      |       |       |       |       |
      | leaf5 | leaf1 | leaf7 | leafx |
      |       |       |       |       |
      +-------+-------+-------+-------+
      |       |       |       |       |
      | leaf3 | brain | feet  | leafy |
      |       |       |       |       |
      +-------+-------+-------+-------+
      |       |       |       |       |
      | leaf8 | leaf2 | leaf6 | leafz |
      |       |       |       |       |
      +-------+-------+-------+-------+

The brain looks in all 8 directions (NW, N, NE, E, SE, S, SW, W) and moves in the direction that reports the smallest distance to a 'CELL'. If none of the 8 directions report a "hit" then use the "hunt" routine to move in a predefined pattern.

The 'feet' check their message buffer and move in the direction indicated. The feet keep moving in the direction indicated until OMOVE fails, then it picks a new direction.

main:

    ; Try to grow and `leaf` in all directions

    ; initialize 1st direction for hunt mode.

    ; (R7, R8) is the most recent successful look vector

brain:

    ; Look in all 8 directions and move toward the one that has the closest thing. If none of them report a "hit" then call the 'hunt' routine to move in a predefined pattern.

direction_array:

    ; Index this array to translate a number (1 thru 8) into an (x,y) vector that can be used with OMOVE.

leaf:

    ; EAT routine. Try to maximize the number of EAT instructions versus the non-eating instructions.

feet:

    ; The feet of this organism continiously check their message buffer for a non-zero value. If it recieves a non-zero message it translates that into a direction and stores the value in (R7,R8)

hunt:

    When all the 8 LOOK operations return no hits, then we must move along some kind of hunting pattern.

    First, If (R7, R8) is non-zero we move in that direction, otherwise we move in the direction specified by R5. If the move using R5 fails we incrment R5 to the next direction.

    Registers used:
     R5 - contains the current direction index
     (R7, R8) most recently sucessful vector

### gonad

This creature is a 13 cell organism WITH HEALING ABILITY. the 
"eat" cells just eat like crazy, while the brain cell LOOK's and
the feet cell MOVE's. The gonad cell tries to reproduce

      +-------+-------+-------+-------+
      |       |       |       |       |
      | eat7  | eat8  | eat9  | eat10 |
      |       |       |       |       |
      +-------+-------+-------+-------+-------+  +--------+
      |       |       |       |       |       |  | new    |
      | eat6  | brain | feet  | eat1  | gonad |  |organism|
      |       |       |       |       |       |  | here   |
      +-------+-------+-------+-------+-------+  +--------+
      |       |       |       |       |
      | eat5  | eat4  | eat3  | eat2  |
      |       |       |       |       |
      +-------+-------+-------+-------+

HEALING:
The healing ability is accomplished by the eat cells. Each eat cell is responsible for growing (if needed) the next eat cell. So eat1 will regrow eat2 if needed.

The order in which this creature is created is as follows:
brain, feet, eat1, eat2, eat3, eat4, eat5, eat6, eat7, eat8, eat9, eat10, gonad

REPRODUCTION:
The gonad cell must reproduce in two steps (two spore creations). If the feet are moving during this time, then there is the possibility of moving between the first and second MAKE-SPORE. To solve this, the MOOD mechanisms will be used to make sure that the feet are not going to be moving for a while.

main:

    ; Look in all 8 directions and move toward the one that has the closest thing. If none of them report a "hit" then call the 'hunt' routine to move in a predefined pattern.

feet:

    ; Feet move the organism. The MOOD of this cell will normally be 0. When the mood is set to 1, this means the gonad's can reproduce with little fear of the organism moving.

gonad:

    ; wait until feet's mood is zero.

    ; wait until feet's mood is non-zero

    ; now we should have enough time to reproduce

direction_array:

    ; Index this array to translate a number (1 thru 8) into an (x,y) vector that can be used with OMOVE.

hunt:

    ; The hunt routine is called by the feet to implement the organism's movement logic.
    When all the 8 LOOK operations return no hits, then we must move along some kind of hunting pattern.

    First, If (R7, R8) is non-zero we move in that direction, otherwise we move in the direction specified by R5. If the move using R5 fails we incrment R5 to the next direction.

    Registers used:
     R5 - contains the current direction index
     (R7, R8) most recently sucessful vector

### greedy_seed

Suggested Environment:

This seed favors very large environments and copious amounts of starting energy. Interesting barriers have no significant effect.

BEHAVIOR:

Register R4 stores the organism's "greediness factor".

This organism will reproduce at the diagonals if its energy is greater than the Greediness.
If its energy goes <= Greediness, it will turn into an attacker that favors organic matter.

If it fails to drop any spores at its current location it will move.

Exercise:  
Make strains with different greediness and set them against one another.

main:

    ; greediness = 6

attacker

sexponddiag:

    ; Place a spore at the diagonals, clockwise.
    ; Store the success flags in registers.
    ; Do not fertilize any of them yet, so that hungry children will not eat this parent.
    ; Fertilize each in turn
    ; Add all the success flags. 
    ; A sum of zero means no reproduction took place.

trymove

Noop

### hanoi

a Tower of Hanoi creature

(c) 2007, Ken Stauffer

Add this creature to a blank universe (no barriers) It will create a bunch of disks and then move them from one side of the screen to the other side.

NOTE: These routines refer to three piles for storing disk.
Piles are encoded as follows:

    -1 =    Left pile
     0 =    Middle pile
     1 =    Right pile

main:

    ; <=== number of disks to play with
    measure_universe
    make_disks
    play_towers_of_hanoi
    ?loop

measure_universe:

    ( -- width height)
    Measure the universe, return the width and height.
    Assumes:
    * Universe is empty, except for itself.
    * No "oval barrier" was used to create the universe.

make_disks:

    ( disks -- )
    Create the initial pile of disk on left-hand side of the universe.

make_disk:

    ( size -- )
    Make a single disk. 'size' is how big the disk is.

put_disk:

    (pile size -- )
    Put a disk down on 'pile'.
    'pile' is where to put the disk.
    'size' is the size of the disk we are putting.

take_disk:

    (pile -- size )
    Pick up a disk from 'pile'
    'pile' is where we will pick up a disk
    'size' is how big of a disk we picked up.
    (pile -- )
    Go to 'pile'.

goto_pile:

    (pile -- )
    Go to 'pile'.

move_disk:

    ( from-pile to-pile -- )
    Move whatever disk is on top of 'from-pile' and place it on top of 'to-pile'.

play_towers_of_hanoi:

    (n src aux dst -- )
    Solve Tower Hanoi problem.
    Implements this algorithm:

    Solve(N, Src, Aux, Dst)
    {
        if N is 0 exit
        Solve(N-1, Src, Dst, Aux)
        Move from Src to Dst
        Solve(N-1, Aux, Src, Dst)
    }

The first invocation of this routine should be:

    N -1 0 1 play_towers_of_hanoi call

(where N is the number of disks)

### hanoi15

Same as `hanoi`

### iamthestar

user contributed creature (doesn't work anymore, but will create a pretty organism, then it just sits there).

main:

    IamSmall
    NewCell
    LiveLife

IamSmall:

    GrowSmall
    SendEnergy

SendEnergy

GrowSmall:

    EatForEver

EatForEver

BeBone

Sex1

Sex3

Sex5

Sex7

LiveLife

NewCell

Expand

Blue_read

Blueprint

### iamthestar2

Very similar to `iamthestar`

### knut

No docu.

### knut2

No docu.

### knut3

No docu.

### seed

a good seed organism

BEHAVIOR:
This creature forever moves in various directions until its forward movement is blocked. As it moves it eats. Before moving in a new direction it will try to reproduce.

main:
    go SOUTH-EAST until blocked (eat along the way)
    go WEST until blocked (eat along the way)
    go SOUTH-WEST until blocked (eat along the way)
    go EAST until blocked (eat along the way)
    go NORTH-EAST until blocked (eat along the way)
    go SOUTH until blocked (eat along the way)
    go NORTH until blocked (eat along the way)
    go NORTH-WEST until blocked (eat along the way)
    keep reproducing between moves
    do it all over again

reproduce:
    compute 1/4 of our energy, and store in R0
    Make 1st spore to the square on our left.
    Put 2nd spore at the same spot (fertilizing it) get the hell out of the way so we don't eat our own babies, or they don't eat us.

### sex_seed

This creature forever moves in various directions until its forward movement is blocked. Before moving in a new direction it will try to reproduce and eat.

The "reproduce" algorithm is both SEXUAL and ASEXUAL.
If there are no organisms nearby, then we reproduce asexually. Otherwise, we try to reproduce sexually (by fertilizing spores, or just creating spores for others to ferilize).

A hybrid sexual/asexual algorithm.

1. If we are touching a spore, then fertilize it
2. Otherwise look around and see if any other organisms can be found.
3. If other organisms are seend, then create a spore (without fertilizing it)
4. Otherwise, create a new organism asexually.

any spores directly touching us?
fertilize any spore we are directly touching
YES, create a lone spore
lay down spore in opposite direction from nearest anything
NO, create organism asexually
move in opposite direction from nearest thing.

### shell_diva

Diva (with protective shell)
Same as 'diva.kf' except is surrounded by a shell of eater leafe
cells, whose sole job is to EAT

      +-------+-------+-------+
      |       |       |       |
      | leaf5 | leaf1 | leaf7 |
      |       |       |       |
      +-------+-------+-------+
      |       |       |       |
      | leaf3 | brain | leaf4 |
      |       |       |       |
      +-------+-------+-------+
      |       |       |       |
      | leaf8 | leaf2 | leaf6 |
      |       |       |       |
      +-------+-------+-------+

Basically `duck_killer.leaf` + `diva`

### sperm

Use this as a seed organism for really interesting evolutionary runs. This organism is really hard to get "started out" But once it gets going, *hold on to the bars*.
This organism requires partners to reproduce. So your initial environment should contain all 8 strains of seed_showoff.kf with equal numbers amounts of energy per strain. Only one will come to dominate later.

main:
    Go to nearest food (eat along the way).
    This code will be repeated for each change in direction.
    go SOUTH until blocked (eat along the way)
    go SOUTH-EAST until blocked (eat along the way)
    go WEST until blocked (eat along the way)
    go SOUTH-WEST until blocked (eat along the way)
    go EAST until blocked (eat along the way)
    go NORTH-EAST until blocked (eat along the way)
    go NORTH until blocked (eat along the way)
    go NORTH-WEST until blocked (eat along the way)
    do it all over again

show_off:
    Organisms will reproduce more if they have an "audience"
    Look around in all 8 directions.
    Reproduce for the number of times you see a *NEARBY* organism.
    Minimum audience distance.

reproduce:
    Stop organisms from reproducing near walls.    
    Minimum wall distance
    compute 1/8 of our energy, and store in R0
    Calculate the escape direction
    Make 1st spore to the square opposite the escape direction
    Put 2nd spore at the same spot (fertilizing it)
    Move in the escape direction

Noop

### sqr

No docu.

### sqr2

No docu.

### sweeper

No docu.

### weasel

(Not working, cell dies at step 40+ after running main. I added `1 ?loop` to keep it alive, but nothing else happens.)

A KFORTH implementation of the "Weasel" program, "The Blind Watchmaker" (by R. Dawkins)

If you thought that program was lacking key features of Darwinism, then you'll even think this program sucks balls (times a thousand) :-)

```
+---+---+---+
| 0 | 1 | 2 |
+---+---+---+
| 3 | X |
+---+---+
```

creature looks at spots 0, 1, 2 and 3 for a spore.

What it sees detmines what letter to draw.

If it doesn't see any spores, then it is the first organism.

main:

    Look in NW-W-SW-S
    letter_draw_routine

letter_draw_routine:
    { monkey_zero }        ; 0000        the first organism
    { draw_m }        ; 0001
    { draw_e }        ; 0010
    { draw_t }        ; 0011
    { draw_h }        ; 0100
    { draw_i }        ; 0101
    { draw_n }        ; 0110
    { draw_k }        ; 0111
    { draw_s }        ; 1000
    { draw_l }        ; 1001
    { draw_a }        ; 1010
    { draw_w }        ; 1011

spore_placement:

    not used
    M
    E

monkey_zero:

    "patient zero", the first organism
    created other organisms, to do the actual drawing.

make_M_monkey
make_E_monkey
make_T_monkey
make_H_monkey
make_I_monkey
make_N_monkey
make_K_monkey
make_S_monkey
make_L_monkey_l
make_A_monkey
make_W_monkey
draw_m
draw_e
draw_t
draw_h
draw_i
draw_n
draw_k
draw_s
draw_l
draw_a
draw_w

### worm

No docu.

### xy_seed

This is a PURE sexual creature. Meaning it doesn't try asexual reproduction in its reproduction routine. It can only be bootstrapped with a small population of other creatures.

(it is also intended to use the "xy " hack)

This creature forever moves in various directions until its forward movement is blocked. Before moving in a new direction it will try to reproduce and eat.

The "reproduce" algorithm is both SEXUAL and ASEXUAL.
If there are no organisms nearby, then we reproduce asexually. Otherwise. we try to reproduce sexually (by fertilizing spores, or just creating spores for others to ferilize).

A pure sexual algorithm.
If we are touching a spore, fertilize it.
Otherwise lay down a spore
reproduce:
    fertilize any spore we are directly touching
    lay down spore in opposite direction from nearest anything
    move in opposite direction from nearest thing.


