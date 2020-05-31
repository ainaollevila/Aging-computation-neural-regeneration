#ifndef Structs_h
#define Structs_h

typedef struct{
    int input,hid1,hid2,output;
}feasible_path;

typedef struct{
    int *inhid1, *hid1hid2, *hid2out;
    int hid1, hid2;
    double *theta_hid1, *theta_hid2, *theta_output; 
    double regeneration_rate;
    
    feasible_path *vec_feas_paths;
    int protected_path;
    int truly_feas_hid1,truly_feas_hid2;
    int *vec_aux_truly_feas_hid1;
    int *vec_aux_truly_feas_hid2;
    int feas_paths;
    int inputs_nonfeas;
} network;

#endif
