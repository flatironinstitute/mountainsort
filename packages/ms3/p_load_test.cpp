/*
 * Copyright 2016-2017 Flatiron Institute, Simons Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "p_load_test.h"

#include <QJsonDocument>
#include <QJsonObject>
#include <QTime>
#include <QDebug>
#include "textfile.h"

bool p_load_test(QString stats_out, P_load_test_opts opts)
{
    QJsonObject stats;
    if (opts.num_cpu_ops) {
        QTime timer;
        timer.start();
        double x = 1;
        for (bigint i = 0; i < opts.num_cpu_ops; i++) {
            x = x + i;
        }
        printf("Print result to force actual computation: %g\n", x);
        stats["cpu_ops_elapsed_msec"] = timer.elapsed();
        double flops = opts.num_cpu_ops / (timer.elapsed() * 1.0 / 1000);
        stats["cpu_megaflops"] = flops / 1e6;
    }
    if (opts.num_read_bytes) {
        QTime timer;
        timer.start();
        FILE* f = fopen("/dev/zero", "rb");
        if (!f) {
            qWarning() << "Unable to open file for reading: /dev/zero";
            return false;
        }
        bigint bufsize = 10000;
        char buffer[bufsize];
        bigint i = 0;
        while (i < opts.num_read_bytes) {
            bigint numread = fread(buffer, 1, bufsize, f);
            (void)numread;
            i += bufsize;
        }
        fclose(f);
        stats["read_elapsed_msec"] = timer.elapsed();
        double mbps = (opts.num_read_bytes * 1.0 / 1e6) / (timer.elapsed() * 1.0 / 1000);
        stats["read_mbps"] = mbps;
    }
    if (opts.num_write_bytes) {
        QTime timer;
        timer.start();
        FILE* f = fopen("/dev/null", "wb");
        if (!f) {
            qWarning() << "Unable to open file for writing: /dev/null";
            return false;
        }
        bigint bufsize = 10000;
        char buffer[bufsize];
        bigint i = 0;
        while (i < opts.num_write_bytes) {
            for (bigint j = 0; j < bufsize; j++)
                buffer[j] = j % 256;
            fwrite(buffer, 1, bufsize, f);
            i += bufsize;
        }
        fclose(f);
        stats["write_elapsed_msec"] = timer.elapsed();
        double mbps = (opts.num_write_bytes * 1.0 / 1e6) / (timer.elapsed() * 1.0 / 1000);
        stats["write_mbps"] = mbps;
    }

    QString json = QJsonDocument(stats).toJson(QJsonDocument::Indented);
    printf("%s\n", json.toUtf8().data());
    return TextFile::write(stats_out, json);
}

