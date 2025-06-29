/******************************************************************************\
 snake.c  Solve the snake:
    "The snake" is 27 (wooden) cubes, each of them except 2 ends
    connected to 2 neighboring cubes, the two ends to only one.
    Each connection is face to face, and can each be rotated.
    The result is that there are straight-line sequences joined
    at right angles in 3 dimensions, each "corner cube" belonging
    to the sequence both before and after it.  The specific lengths
    of these straight-line sequences, in order from one end to the
    other, and thereby double-counting the corner cubes, is:
           3, 3, 3, 3, 2, 2, 2, 3, 3, 2, 2, 3, 2, 3, 2, 2, 3

    The goal is to manipulate the snake into a 3x3 "big cube" while
    preserving these constraints.

    This program uses a brute force approach. Whenever we establish
    a location for one of those segments, we'll need to look in all 4
    perpendicular directions from its end, and not only try to put the
    next segment there (if there's room, and no other segment obstructs
    it), but also do the same at the end of THAT segment, etc. If/when
    we can't move forward this way, we need to backtrack to the last
    decision point - that is, the last segment where there's a
    perpendicular direction we haven't tried, and then try that. We can
    envision all these choices as a big decision tree, and can use a
    stack - specifically, the function call stack - to keep track of
    where we are, moving away from the root (to try more segments and
    directions) by calling more functions, or backtracking back towards
        the root (when failing to proceed) by returning from those functions.

    We'll cheat a bit at the beginning. Since the first 4 segments in
    the snake are 3 long, we know that the first segment must end up
    along a corner of the "big cube". So instead of starting out by
    trying all kinds of directions for the first segment, we'll just
    start at a corner of the big cube and try (and succeed) to lay the
    first segment down from there.

    A note about "directions": A direction will designate both an axis in
    the big (3x3x3 array) cube (x, y, or z, which correspond roughly
    to dimensions 0, 1, and 2) and whether we're moving in a
    negative or positive direction along that axis. It is easiest for
    us to consider this as a unit vector in 3D, that is, a 3-element 1D
    array (called offset) where all the elements are 0 except one which
    is either a +1 or -1. We'll use a separate variable, dim, to designate
    which of the elements of offset is non-zero.

    So, primary data structures are:
        int cube[3][3][3] represents the big cube, non-zero element
             means it's full (we'll put (segment # + 1) in there).
        int offset[3] and int dim represent a direction (as above).
        int snake[] is the sequence of segment lengths in the snake.
        int pos[3] is coordinate of a position in cube.
        int segment is the index into snake of segment we're working on.

    The two principal functions we'll use are:
        int go_dir(dim, offset, cube, pos, segment)
            Tries to put segment into cube starting at
              pos in direction offset (in dimension dim).
            Returns a 1 if it succeeds (ie. there's room, it's open,
              and all following segments also fit) or a 0 if not.

        int find_new_dir(cube, segment, pos, not_dim)
            Tries every direction (which is not along dimension
               not_dim) in cube for segment, starting at pos. 
            Returns a 1 if it finds one (that works all the way
               from here on out) or a 0 if not.
\******************************************************************************/

#include <stdio.h>
#include <string.h>

int go_dir();
int find_new_dir();

/* OFFSET is a simple macro to give the location in cube x which is n "offset"
    units from pos */

#define    OFFSET(x,n)    x[pos[0]+offset[0]*n]        \
             [pos[1]+offset[1]*n]        \
             [pos[2]+offset[2]*n]

/* Make WATCH_BACKTRACKS non-zero to see it trying backtracks */
#define WATCH_BACKTRACKS    0

/* Segment lengths in the snake */

static int snake[]= {3, 3, 3, 3, 2, 2, 2, 3, 3, 2, 2, 3, 2, 3, 2, 2, 3};


int main() {
    /* initial 3-dim cube, start with pos filled with beg of segment 0 */

    int cube[3][3][3] = { {{1, 0, 0},{0, 0, 0}, {0, 0, 0}},
                  {{0, 0, 0},{0, 0, 0}, {0, 0, 0}},
                  {{0, 0, 0},{0, 0, 0}, {0, 0, 0}}};

    /* Start with segment 0, pos 0 0 0 in with vector offset 1, 0, 0. */

    int pos[3] = {0, 0, 0};
    int offset[3] = {1, 0, 0};
    int segment = 0;

    /* Lay down first segment, and all after */

    if (!go_dir(0, offset, cube, pos, 0)) {
    
        /* Failed overall ?! Shouldn't happen, snake too long? */

        printf("Couldn't do it, not enough room!\n");
    }
}

/* Find a new direction (starting at position pos, for segment, in some non-dim
 * dimension) that works from there on out. */

int find_new_dir(cube, segment, pos, not_dim)
int cube[3][3][3];
int segment;
int pos[3];
int not_dim;
{
    int offset[3] = {0, 0, 0};
    int dim;

    /* If this segment is past the end of the snake, we made it! */

    if (segment >= (sizeof(snake) / sizeof(int))) {

    /* Print out cube, containing segment numbers */

    printf("Final cube, 3 sequential slices, segment # filling each loc\n");
    for (int i=0; i <=2; i++) {

        printf("Slice %d\n", i);

        for (int j=0; j <=2; j++) {
            for (int k=0; k <=2; k++) {
            printf(" %2d", cube[i][j][k]);
            }
            printf("\n");
        }
                printf("\n");
    }

    /* Prepare for directions to spill out, then return with success */

    printf("Segment directions in reverse order:\n");

    return(1);
    }

    /* OK, new segment. For every dimension not equal to notdim... */

    for (dim = 0; dim <= 2; dim++) {
    if (dim != not_dim) {

        /* If positive in that dimension works, we made it! */

        offset[dim] = 1;
        if (go_dir(dim, offset, cube, pos, segment)) return (1);

        /* If negative in that dimension works, we made it! */

        offset[dim] = -1;
        if (go_dir(dim, offset, cube, pos, segment)) return (1);

        /* Nope, nothing in this dimension */

        offset[dim] = 0;
    }
    }

    /* Nope, nothing in any dimension, fail. */

    return(0);
}

/* Try to go segment distance in offset direction from position pos */

int go_dir(dim, offset, cube, pos, segment)
int dim;
int offset[3];
int cube[3][3][3];
int pos[3];
int segment;
{
    int newcube[3][3][3];
    int newpos[3];

    int distance = snake[segment] - 1;

#if WATCH_BACKTRACKS
    printf("going to try segment %d length %d from %d %d %d dir %d\n",
    segment, distance, pos[0],pos[1],pos[2],(dim+1)*offset[dim]);
#endif

    /* Fail if not enough room from pos to edge in direction */

    if (offset[dim] < 0) {
    if (pos[dim] - distance < 0) return(0);
    } else {
    if (pos[dim] + distance > 2) return(0);
    }
  
    /* Fail if any locations in path are already full */

    if (OFFSET(cube,1) || (distance == 2 && OFFSET(cube,2))) return(0);
         
    /* We're committing to try this direction. Make new cube with segment in */

    memcpy(newcube, cube, sizeof(int)*3*3*3);
    OFFSET(newcube,1) = segment+1;
    if (distance == 2) OFFSET(newcube,2) = segment+1;

    /* Move to end */

    memcpy(newpos, pos, sizeof(int)*3);
    newpos[dim] += distance * offset[dim];

    /* Still, fail if we can't find a direction from there that works */

    if (!find_new_dir(newcube, segment+1, newpos, dim)) return(0);

    /* Yay, it worked! Report what we did on this level, and return  success */

    printf("Segment of %d in dir %d from position %d %d %d to %d %d %d\n",
         distance + 1, (dim + 1) * offset[dim], pos[0], pos[1], pos[2],
                newpos[0],newpos[1],newpos[2]);
    return(1);
}
 