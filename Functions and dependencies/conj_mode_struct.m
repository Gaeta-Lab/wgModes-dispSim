function conjMode=conj_mode_struct(mode)

mode_E=mode.E;
mode_H=mode.H;

conjMode.E=conj(mode_E);
conjMode.H=conj(mode_H);