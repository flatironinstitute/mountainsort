#include "p_create_firings.h"

#include "mda.h"

#include <mda32.h>

bool p_create_firings(QString event_times, QString labels, QString amplitudes, QString firings_out, int central_channel)
{
    Mda ET(event_times);
    Mda LA(labels);
    Mda32 AM;

    bigint L = ET.totalSize();

    int R = 3;
    if (!amplitudes.isEmpty()) {
        R = 4;
        AM.read(amplitudes);
    }
    Mda firings(R, L);
    for (bigint i = 0; i < L; i++) {
        firings.setValue(central_channel, 0, i);
        firings.setValue(ET.value(i), 1, i);
        firings.setValue(LA.value(i), 2, i);
        if (R >= 4)
            firings.setValue(AM.value(i), 3, i);
    }

    return firings.write64(firings_out);
}
