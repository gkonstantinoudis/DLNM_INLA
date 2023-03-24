# DLNM_INLA
Temperature and mortality is a U-shape. We focus on summer mortality and want to fit the following models with R-INLA:
- Linear threshold model in which the slope above the threshold varies in time and space independently
- RW model for temperature that has a spatial and temporal structure (Q_rw x Q_space and Q_rw x Q_time)
- We can expand both approaches to let them have a lag dimention, i.e. if lag: iid then we look at smth like this Q_iid x Q_rw x Q_space and Q_iid x Q_rw x Q_time
