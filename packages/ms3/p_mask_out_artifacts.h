/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef MASK_OUT_ARTIFACTS_H
#define MASK_OUT_ARTIFACTS_H

#include <QString>

bool p_mask_out_artifacts(const QString& timeseries_path, const QString& timeseries_out_path, double threshold, int interval_size);

#endif // MASK_OUT_ARTIFACTS_H
