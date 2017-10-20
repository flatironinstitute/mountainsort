#ifndef P_LOAD_TEST_H
#define P_LOAD_TEST_H

#include "mlutil.h"

struct P_load_test_opts {
    bigint num_cpu_ops = 0;
    bigint num_read_bytes = 0;
    bigint num_write_bytes = 0;
};

bool p_load_test(QString stats_out, P_load_test_opts opts);

bool p_misc_test(QString dir, QString info_out);

#endif // P_LOAD_TEST_H
