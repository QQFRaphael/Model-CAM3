#include "ipc.h"

program SCAM_PROG {
    version SCAM_VERS {
        int               INIT_MODEL(struct init_data)               = 1;
        field_info        GET_FIELD_INFO(int)                        = 2;
        field_data        GET_FIELD(int)                             = 3;
        level_data        GET_LEVELS(void)                           = 4;
	dataset_var_real  GET_DATASET_VAR_REAL(struct dataset_req)   = 5;
	dataset_var_int   GET_DATASET_VAR_INT(struct dataset_req)    = 6;
        int               GET_DATASET_DIMSIZE(struct dataset_req)    = 7;
        int               SET_FIELD(struct field_data)               = 8;
	int               RESET_FIELD(int)                           = 9;
        int               STEP(void)                                 = 10;
    } = 1;
} = 0x3ffffff1;








