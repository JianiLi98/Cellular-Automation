/* Program to perform 1D cellular automaton (CA) computations and to use 1D CA
   to solve the density classification problem.

  Skeleton program written by Artem Polyvyanyy, http://polyvyanyy.com/,
  September 2024, with the intention that it be modified by students
  to add functionality, as required by the assignment specification.
  All included code is (c) Copyright University of Melbourne, 2024.

  Authorship Declaration:

  (1) I certify that except for the code provided in the initial skeleton file,
  the program contained in this submission is completely my own individual
  work, except where explicitly noted by further comments that provide details
  otherwise. I understand that work that has been developed by another student,
  or by me in collaboration with other students, or by non-students as a result
  of request, solicitation, or payment, may not be submitted for assessment in
  this subject. I understand that submitting for assessment work developed by
  or in collaboration with other students or non-students constitutes Academic
  Misconduct, and may be penalized by mark deductions, or by other penalties
  determined via the University of Melbourne Academic Honesty Policy, as
  described at https://academicintegrity.unimelb.edu.au.

  (2) I also certify that I have not provided a copy of this work in either
  softcopy or hardcopy or any other form to any other student, and nor will I
  do so until after the marks are released. I understand that providing my work
  to other students, regardless of my intention or any undertakings made to me
  by that other student, is also Academic Misconduct.

  (3) I further understand that providing a copy of the assignment specification
  to any form of code authoring or assignment tutoring service, or drawing the
  attention of others to such services and code that may have been made
  available via such a service, may be regarded as Student General Misconduct
  (interfering with the teaching activities of the University and/or inciting
  others to commit Academic Misconduct). I understand that an allegation of
  Student General Misconduct may arise regardless of whether or not I personally
  make use of such solutions or sought benefit from such actions.

  Signed by: Jiani Li 1569549
  Dated:     September 29th 2024
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* #DEFINE'S -----------------------------------------------------------------*/
#define SDELIM "==STAGE %d============================\n"   // stage delimiter
#define MDELIM "-------------------------------------\n"    // delimiter of -'s
#define THEEND "==THE END============================\n"    // end message

#define CRTRNC '\r'     // carriage return character
#define NEW_LINE '\n'   // new line character
#define PRINT_NL "\n"   // for printing new line
#define NBRHDS 8        // number of possible neighborhoods

#define ONE_NBRHD 3     // number of cells per neighborhood
#define INIT_TIME 0     // initial state's time step value

#define ON 1            // denotes state of an on cell
#define OFF 0           // denotes state of an off cell
#define ON_CHAR '*'     // character representation of on cell
#define OFF_CHAR '.'    // character representation of off cell

#define BIN_BASE 1   // the zero-th position's decimal value for binary numbers
#define BIN_FAC 2    // the multiplication factor to the next position

#define PREV_CLL 0   // index position for a previous cell in a neighborhood
#define CURR_CLL 1   // index position for current cell in a neighborhood
#define NEXT_CLL 2   // index position for the next cell in a neighborhood

#define RULE_184 184 // the first update rule and...
#define RULE_232 232 // second update rule for density classification problem

#define ALL_SAME 1   // sequence has all on or all off cells
#define NOT_SAME 0   // sequence has not only one state of cells

#define S_ZERO 0     // Stage 0
#define S_ONE 1      // Stage 1
#define S_TWO 2      // Stage 2

/* TYPE DEFINITIONS ----------------------------------------------------------*/
typedef char cells_t;            // base type to store states of cells
typedef struct state state_t;           // a cellular automaton state
typedef unsigned char rule_t[NBRHDS];   // an elementary CA update rule function
typedef unsigned char nbrhd_t[ONE_NBRHD]; // a neighborhood of cells

struct state {                   // a state in a CA is defined by
    cells_t*        clls;        // ... an array of cells and
    state_t*        next;        // ... a link to the next state
};

typedef struct {                 // a run of a CA consists of
    state_t*        init;        // ... the initial state and
    state_t*        curr;        // ... the current state,
} run_t;                         // implemented as a linked list of states


typedef struct {                 // an elementary CA is defined by
    unsigned char   code;        // ... a code of the update rule,
    unsigned int    size;        // ... a number of cells,
    unsigned int    elts;        // ... number of elements to store cells,
    unsigned int    time;        // ... the current time step,
    rule_t          rule;        // ... an update rule function, and
    run_t*          run;         // ... a run of state steps
} CA_t;

/* USEFUL FUNCTIONS ----------------------------------------------------------*/
int             mygetchar(void);                // getchar() that skips
                                                //    carriage returns
CA_t* make_new_CA(void);
void do_stage_0(CA_t* CA);
void get_bin_array(unsigned char* array, int decimal, int limit);
void print_rule_header(void);
void print_rule(rule_t rule, int limit);
void print_state(state_t *state, int time_step, int size);
run_t* make_empty_run(void);
state_t* get_init_state(int size);
CA_t* initialise_automation(CA_t* CA, state_t* init_state);
void do_stage_1(CA_t* CA);
void get_neighbors(CA_t* CA, int pos, nbrhd_t neighborhood);
cells_t evolve(nbrhd_t neighborhood, int limit, rule_t rule);
CA_t* evolve_automation(CA_t* CA);
CA_t* run_automation(CA_t* CA, int start, int end);
state_t* get_state(CA_t* CA, int target_time);
int  get_on_num(CA_t* CA, int cll_pos, int start_time);
void count_cll_states(CA_t* CA, unsigned int tot_time);
int  get_time_steps(unsigned int size, int rule);
void do_stage_2(CA_t* CA, int n, int m);
void print_S2_rules(int rule, int times);
void answer_density_prob(CA_t* CA, unsigned int times_num);
int  are_unequal_clls(state_t* state, int size);
void free_automation(CA_t* CA);
void print_delimiter(void);

/* WHERE IT ALL HAPPENS ------------------------------------------------------*/
int main(int argc, char *argv[]) {
    // Message from Artem: The proposed in this skeleton file #define's,
    // typedef's, and struct's are the subsets of those from my sample solution
    // to this assignment. You can decide to use them in your program, or if
    // you find them confusing, you can remove them and implement your solution
    // from scratch. I will share my sample solution with you at the end of
    // the subject.

    /* create and initialise the cellular automation */
    CA_t* automata = make_new_CA();
    scanf("%u\n%hhu\n", &(automata->size), &(automata->code));

    printf(SDELIM, S_ZERO);
    do_stage_0(automata);

    scanf("%u\n", &(automata->time)); 
    printf(SDELIM, S_ONE);  
    do_stage_1(automata);

    int n = get_time_steps(automata->size, RULE_184);
    int m = get_time_steps(automata->size, RULE_232);
    printf(SDELIM, S_TWO);
    do_stage_2(automata, n, m);

    printf(THEEND);

    /* free the automation */
    free_automation(automata);

    return EXIT_SUCCESS;        // algorithms are fun!!!
}

/* USEFUL FUNCTIONS ----------------------------------------------------------*/

/* Makes a new cellular automaton and returns a pointer to it.
*/
CA_t* make_new_CA(void) {
    CA_t* new_CA = (CA_t *)malloc(sizeof(*new_CA));
    assert(new_CA);
    // initialise the details of the automation
    new_CA->code = new_CA->size = new_CA->time = 0;
    return new_CA;
}

/* Performs Stage 0 with tasks split between functions.
*/
void do_stage_0(CA_t* CA) {
    printf("SIZE: %d\n", CA->size);
    printf("RULE: %d\n", CA->code);
    print_delimiter();

    /* print the rule header and the rule in required format */
    print_rule_header();
    get_bin_array(CA->rule, CA->code, NBRHDS);
    print_rule(CA->rule, NBRHDS);
    print_delimiter();
    
    /* initialiise the CA with the initial state specified from stdin */
    CA->run = make_empty_run();
    state_t* init_state = get_init_state(CA->size);
    CA = initialise_automation(CA, init_state);
    print_state(CA->run->init, INIT_TIME, CA->size);
}

/* Takes a decimal integer and converts it into a binary number with each 
   digit stored in a position in the array.
*/
void get_bin_array(unsigned char* array, int decimal, int limit) {
    int position = 0;
    /* fill the array until the limit provided */
    while(position < limit) {
        if(decimal == 0) {
            /* fill the rest with 0s if decimal is done being divided by 2 */
            array[position ++] = 0;
        } else {
            /* keep filling in the remainder of the decimal divided by 2 */
            array[position ++] = decimal % 2;
            decimal = decimal / 2;
        }
    }
}

/* Helps the visualisation of the update rule by printing out a rule header 
   of all the possible neighborhoods in order.
*/
void print_rule_header(void) {
    nbrhd_t neighborhood;
    for(int nbrhd = 0; nbrhd < NBRHDS; nbrhd ++) {
        /* gets the binary array of the neighborhood */
        get_bin_array(neighborhood, nbrhd, ONE_NBRHD);
        printf(" ");
        for(int pos = ONE_NBRHD-1; pos >= 0; pos --) {
            printf("%d", neighborhood[pos]);
        }
    }
    printf(PRINT_NL);
}

/* Prints the rule aligned with the header for visualisation. 
*/
void print_rule(rule_t rule, int limit) {
    for(int pos = 0; pos < limit; pos ++) {
        printf("%3d ", rule[pos]);
    }
    printf(PRINT_NL);
}

/* This function is an adapted version of the make_empty_list 
   function by Alistair Moffat: https://people.eng.unimelb.edu.au/
   ammoffat/ppsaa/c/listops.c
   Makes an empty run of cellular automation.
*/
run_t* make_empty_run(void) {
    run_t* run = (run_t *)malloc(sizeof(*run));
    assert(run);
    run->curr = run->init = NULL;
    return run;
}

/* Reads from stdin the specified initial state of the cellular automation
   and returns a pointer to the state, provided the size of cells.
*/
state_t* get_init_state(int size) {
    int ch, clls_num = 0;

    state_t* init_state = (state_t *)malloc(sizeof(*init_state));
    assert(init_state);
    init_state->clls = (cells_t *)malloc(size * sizeof(cells_t));
    assert(init_state->clls);
    init_state->next = NULL;

    /* while have not reached a new-line character */
    while((ch=mygetchar()) != NEW_LINE) {
        /* stores the input initial state denoted by characters */
        init_state->clls[clls_num++] = ch;
    }
    return init_state;
}

/* Helps to initialise a run of the automation with an initial state provided.
*/
CA_t* initialise_automation(CA_t* CA, state_t* init_state) {
    CA->run->init = CA->run->curr = init_state;
    return CA;
}

/* Performs Stage 1: evolve the automation for the number of time steps 
   specified in stdin with helps from other functions.
*/
void do_stage_1(CA_t* CA) {
    print_state(CA->run->init, INIT_TIME, CA->size);
    /* run the CA for the specified times and print the results */
    run_automation(CA, INIT_TIME+1, CA->time);
    /* prints the requested cell's counts of on and off states */
    count_cll_states(CA, CA->time);
}

/* Takes the CA and the requested cell position and returns the states of 
   the neighborhood. Allows edge cases which wraps around the array. */
void get_neighbors(CA_t* CA, int pos, nbrhd_t neighborhood) {
    int size = CA->size;
    /* assign the state of each cell and ensure to handle edge cells */
    neighborhood[PREV_CLL] = CA->run->curr->clls[(pos - 1 + size) % size];
    neighborhood[CURR_CLL] = CA->run->curr->clls[pos];
    neighborhood[NEXT_CLL] = CA->run->curr->clls[(pos + 1) % size];
}

/* Takes one neighborhood denoted in characters, evaluates its evolution
   based on the update rule, and returns as a character. */
cells_t evolve(nbrhd_t neighborhood, int limit, rule_t rule) {
    int decimal = 0, result = 0, pos;
    int base_factor = BIN_BASE;
    
    /* starting from the right-most character of the neighborhood */
    for(pos = limit-1; pos >= 0; pos --) {
        int curr = ON;
        /* determine the decimal representation of the cell */
        if(neighborhood[pos] == OFF_CHAR) {
            curr = OFF;
        }
        /* convert neighborhood of cells to one decimal integer */
        decimal += curr * base_factor;
        base_factor *= BIN_FAC;
    }
    
    /* use the decimal to index the correct evolved result */
    result = rule[decimal];
    if(result == 1) {
        return ON_CHAR;
    } else {
        return OFF_CHAR;
    }
}

/* Evolve the cellular automation once by adding a new state to the run. 
*/
CA_t* evolve_automation(CA_t* CA) {
    state_t* new = (state_t *)malloc(sizeof(*new));
    assert(new);
    new->clls = (cells_t *)malloc(CA->size * sizeof(cells_t));
    assert(new->clls);
    new->next = NULL;
    
    /* iterate through the cells of the current state */
    for(int pos = 0; pos < CA->size; pos ++) {
        /* get the current neighborhood to be evolved */
        nbrhd_t neighborhood;
        get_neighbors(CA, pos, neighborhood);
        /* add the evolved cell to the new state */
        new->clls[pos] = evolve(neighborhood, ONE_NBRHD, CA->rule);
    }

    /* add the new state to the run */
    CA->run->curr->next = new;
    CA->run->curr = new;
    
    return CA;
}

/* Takes a start and an end time step, performs a run of the cellular 
   automation based on the update rule. Prints each evolution. 
*/
CA_t* run_automation(CA_t* CA, int start, int end) {
    /* ensures a correct update rule based on the code */
    get_bin_array(CA->rule, CA->code, NBRHDS);

    /* evolve for the requested times */
    for(int time = start; time < end+1; time ++) {
        CA = evolve_automation(CA);
        print_state(CA->run->curr, time, CA->size);
    }
    print_delimiter();
    
    return CA;
}

/* Takes a requested time step and returns a pointer to the 
   corresponding state. 
*/
state_t* get_state(CA_t* CA, int target_time) {
    state_t* target_state = CA->run->init;
    
    /* go through the run until reach the requested time step */
    for(int time = 0; time < target_time; time ++) {
        target_state = target_state->next;
    }
    
    return target_state;
}

/* Takes the requested cell position and the start time, returns a 
   count of the on states of the run. 
*/
int get_on_num(CA_t* CA, int cll_pos, int start_time) {
    int on_num = 0;
    state_t* target_state = get_state(CA, start_time);

    /* iterate through each state until the end of the run */
    while(target_state != NULL) {
        /* increment the count of on states when encountered */
        if(target_state->clls[cll_pos] == ON_CHAR) {
            on_num ++;
        }
        target_state = target_state->next;
    }

    return on_num;
}

/* Performs and prints the counts of on and off states of requested cell 
   based on inputs of instructions from stdin. 
*/
void count_cll_states(CA_t* CA, unsigned int tot_time) {
    int cll_pos, start_time;
    scanf("%d,%d\n", &cll_pos, &start_time); /* get the instructions */

    int on_num = 0, off_num = 0;
    /* gets the on state count */
    on_num = get_on_num(CA, cll_pos, start_time);
    /* calculate the off state count */
    off_num = tot_time - (start_time - 1) - on_num;
    printf("#ON=%d #OFF=%d CELL#%d START@%d\n", 
        on_num, off_num, cll_pos, start_time);
}

/* Calculates the execution time steps for Fuks' update rules.
*/
int get_time_steps(unsigned int size, int rule) {
    int times;
    if(rule == RULE_184) {
        times = (int)floor((size - 2) / 2);
        return times;  
    } else {
        times = (int)floor((size - 1) / 2);
        return times;
    }
}

/* Performs Stage 2 with tasks split between functions.
*/
void do_stage_2(CA_t* CA, int n, int m) {
    int curr_time = CA->time;
    
    /* evolve for rule 184 */
    print_S2_rules(RULE_184, n);
    print_state(CA->run->curr, curr_time, CA->size);
    CA->code = RULE_184;
    CA = run_automation(CA, curr_time+1, curr_time+n);
    curr_time += n;

    /* evolve for rule 232 */
    print_S2_rules(RULE_232, m);
    print_state(CA->run->curr, curr_time, CA->size);
    CA->code = RULE_232;
    CA = run_automation(CA, curr_time+1, curr_time+m);
    curr_time += m;
    
    /* prints the requested cell's counts of on and off states */
    count_cll_states(CA, curr_time);
    print_delimiter();
    
    /* answer the density classification problem */
    state_t* problem_automation = get_state(CA, CA->time);
    print_state(problem_automation, CA->time, CA->size);
    answer_density_prob(CA, CA->time);
}

/* Prints the rules and time steps for Stage 2. */
void print_S2_rules(int rule, int times) {
    printf("RULE: %d; STEPS: %d.\n", rule, times);
    print_delimiter();
}

/* Determines and prints whether there are more on cells, more off 
   cells, or equal number of on and off cells. 
*/
void answer_density_prob(CA_t* CA, unsigned int times_num) {
    /* if the final pattern is uniform */
    if(are_unequal_clls(CA->run->curr, CA->size)) {
        if(CA->run->curr->clls[0] == ON_CHAR) {
            printf("AT T=%d: #ON/#CELLS > 1/2\n", times_num);
        } else if(CA->run->curr->clls[0] == OFF_CHAR) {
            printf("AT T=%d: #ON/#CELLS < 1/2\n", times_num);
        }
    /* if the final pattern presents an alternating sequence */
    } else {
        printf("AT T=%d: #ON/#CELLS = 1/2\n", times_num);
    }
}

/* Takes a state and size of the automation, returns ALL_SAME if the final
   sequence is consisted of all on or all off cells.
*/
int are_unequal_clls(state_t* state, int size) {
    int all_same = ALL_SAME;
    cells_t cll = state->clls[0];

    /* check if the state is alternating starting from the second cell */
    for(int pos = 1; pos < size; pos ++) {
        if(state->clls[pos] != cll) {
            /* the sequence is alternating */
            all_same = NOT_SAME;
            break;
        }
    }
    return all_same;
}

/* Adapted version of the free_list function by
   Alistair Moffat: https://people.eng.unimelb.edu.au/ammoffat/ppsaa/c/listops.c

   Free a cellular automation: this is adapted to free all the states, arrays
   of cells in the run, then free the run, and finally the automation. 
*/
void free_automation(CA_t* CA) {
    state_t *curr_state, *prev_state;
    /* start from the first state */
    curr_state = CA->run->init;

    /* while it has not reached the end of the run */
    while(curr_state != NULL) {
        prev_state = curr_state;
        curr_state = curr_state->next;
        /* free the array cells in the state first */
        if(prev_state->clls != NULL) {
            free(prev_state->clls);
            prev_state->clls = NULL;
        }
        /* free the state it self */
        free(prev_state);
    }

    /* free the run of the automation */
    if(CA->run != NULL) {
        free(CA->run); 
        CA->run = NULL;
    }
    /* free the whole automation */
    free(CA);
    CA = NULL;
}

/* Prints the provided state and time step in correct format. 
*/
void print_state(state_t *state, int time_step, int size) {
    printf("%4d: ", time_step);
    /* prints the array of cells denoted by characters */
    for(int pos = 0; pos < size; pos ++) {
        printf("%c", state->clls[pos]);
    }
    printf(PRINT_NL);
}

void print_delimiter(void) {
    printf(MDELIM);
}

// An improved version of getchar(); skips carriage return characters.
// NB: Adapted version of the mygetchar() function by Alistair Moffat
int mygetchar() {
    int c;
    while ((c=getchar())==CRTRNC);          // skip carriage return characters
    return c;
}

/* THE END -------------------------------------------------------------------*/