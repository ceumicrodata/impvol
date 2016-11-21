th = 4;

t1 = ['table_th_', num2str(th), '_lac_inf'];

load_model_specifications(t1)

run_model([t1, '_baseline'])

run_model([t1, '_nosectoral'])

run_model([t1, '_kappa1972'])

run_model([t1, '_kappa1972_nosectoral'])

create_table(t1)



t2 = ['table_th_', num2str(th), '_lac_2000'];

load_model_specifications(t2)

run_model([t2, '_baseline'])

run_model([t2, '_nosectoral'])

run_model([t2, '_kappa1972'])

run_model([t2, '_kappa1972_nosectoral'])

create_table(t2)
