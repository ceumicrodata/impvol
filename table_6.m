th = 4;

t1 = ['table_th_', num2str(th), '_lac_inf_bt_china_1'];

load_model_specifications(t1)

run_model([t1, '_baseline'])

run_model([t1, '_nosectoral'])

run_model([t1, '_kappa1972'])

run_model([t1, '_kappa1972_nosectoral'])

create_table(t1)



th = 4;

t1 = ['table_th_', num2str(th), '_lac_inf_bt_china_2'];

load_model_specifications(t1)

run_model([t1, '_baseline'])

run_model([t1, '_nosectoral'])

run_model([t1, '_kappa1972'])

run_model([t1, '_kappa1972_nosectoral'])

create_table(t1)