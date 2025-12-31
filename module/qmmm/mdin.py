from pathlib import Path
from typing import List, Dict, Generator
from functools import lru_cache

class MdinGenerator:
    __slots__ = ('_templates',)

    def __init__(self):
        self._templates = {}
        self._init_templates()

    def _init_templates(self):
        self._templates['min'] = 'Minmize all the hydrogens\n &cntrl\n ifqnt=1,\n imin=1,           ! Minimize the initial structure\n maxcyc=1,          ! Maximum number of cycles for minimization\n ncyc=1,\n ntb=1,            ! Constant volume\n ntp=0,            ! No pressure scaling\n ntf=1,            ! Complete force evaluation\n ntwx= 1,       ! Write to trajectory file every ntwx steps\n ntpr= 1,       ! Print to mdout every ntpr steps\n ntwr= 1,       ! Write a restart file every ntwr steps\n cut=  10.0,        ! Nonbonded cutoff in Angstroms\n ntr=1,            ! Turn on restraints\n restraintmask=":1-314 & !@H=", ! atoms to be restrained\n restraint_wt=1.0, ! force constant for restraint\n ntxo=1,           ! Write coordinate file in ASCII format\n ioutfm=0,         ! Write trajectory file in ASCII format\n /\n &qmmm\n qmmask = \':111,55\',\n qmcharge = 0,\n qm_theory = \'XTB\',\n spin = 2,\n verbosity = 2,\n printcharges = 1,\n qmcut = 10.0,\n qm_ewald = 1,\n writepdb = 1,\n /\n &xtb\n qm_level = \'GFN2-xTB\',\n tfermi = 300.0,\n accuracy = 1.d-4\n maxiter = 250\n mmhardness = 0.0\n debug = F\n /\n'
        self._templates['sp'] = 'SP calculation\n &cntrl\n ifqnt=1,\n imin=6,           ! Minimize the initial structure\n nstlim=1,\n dt=0.0,\n ntb=0,            ! Constant volume\n ntp=0,            ! No pressure scaling\n ntf=1,            ! Complete force evaluation\n ntwx= 1,       ! Write to trajectory file every ntwx steps\n ntpr= 1,       ! Print to mdout every ntpr steps\n ntwr= 1,       ! Write a restart file every ntwr steps\n cut=  10.0,        ! Nonbonded cutoff in Angstroms\n ntr=1,            ! Turn on restraints\n restraintmask=":1-314 & !@H=", ! atoms to be restrained\n restraint_wt=1.0, ! force constant for restraint\n ntxo=1,           ! Write coordinate file in ASCII format\n ioutfm=0,         ! Write trajectory file in ASCII format\n /\n &qmmm\n qmmask = \':111,55\',\n qmcharge = ,\n qm_theory = \'XTB\',\n spin = 2,\n verbosity = 2,\n printcharges = 1,\n qmcut = 10.0,\n qm_ewald = 1,\n writepdb = 1,\n /\n &xtb\n qm_level = \'GFN2-xTB\',\n tfermi = 300.0,\n accuracy = 1.d-4\n maxiter = 250\n mmhardness = 0.0\n debug = F\n /\n'
        self._templates['vdwm'] = 'MMSP calculation\n &cntrl\n imin=6,           ! Minimize the initial structure\n nstlim=1,\n dt=0.0,\n ntb=0,            ! Constant volume\n ntp=0,            ! No pressure scaling\n ntf=1,            ! Complete force evaluation\n ntwx= 1,       ! Write to trajectory file every ntwx steps\n ntpr= 1,       ! Print to mdout every ntpr steps\n ntwr= 1,       ! Write a restart file every ntwr steps\n cut=  10.0,        ! Nonbonded cutoff in Angstroms\n ntr=1,            ! Turn on restraints\n restraintmask=":1-314 & !@H=", ! atoms to be restrained\n restraint_wt=1.0, ! force constant for restraint\n ntxo=1,           ! Write coordinate file in ASCII format\n ioutfm=0,         ! Write trajectory file in ASCII format\n /\n'
        sys2_templates = {'min_sys2': 'Min SYS2 calculation\n     &cntrl\n     ifqnt=1,\n     imin=1,           ! Minimize the initial structure\n     maxcyc=1,          ! Maximum number of cycles for minimization\n     ncyc=1,\n     ntb=1,            ! Constant volume\n     ntp=0,            ! No pressure scaling\n     ntf=1,            ! Complete force evaluation\n     ntwx= 1,       ! Write to trajectory file every ntwx steps\n     ntpr= 1,       ! Print to mdout every ntpr steps\n     ntwr= 1,       ! Write a restart file every ntwr steps\n     cut=  10.0,        ! Nonbonded cutoff in Angstroms\n     ntr=1,            ! Turn on restraints\n     restraintmask=":1-314 & !@H=", ! atoms to be restrained\n     restraint_wt=1.0, ! force constant for restraint\n     ntxo=1,           ! Write coordinate file in ASCII format\n     ioutfm=0,         ! Write trajectory file in ASCII format\n     /\n     &qmmm\n     qmmask = \':55\',\n     qmcharge = 0,\n     qm_theory = \'XTB\',\n     spin = 2,\n     verbosity = 2,\n     printcharges = 1,\n     qmcut = 10.0,\n     qm_ewald = 1,\n     /\n     &xtb\n     qm_level = \'GFN2-xTB\',\n     tfermi = 300.0,\n     accuracy = 1.d-4\n     maxiter = 250\n     mmhardness = 0.0\n     debug = F\n     /\n    ', 'sp_sys2': 'SP SYS2 calculation\n     &cntrl\n     ifqnt=1,\n     imin=6,           ! Minimize the initial structure\n     nstlim=1,\n     dt=0.0,\n     ntb=0,            ! Constant volume\n     ntp=0,            ! No pressure scaling\n     ntf=1,            ! Complete force evaluation\n     ntwx= 1,       ! Write to trajectory file every ntwx steps\n     ntpr= 1,       ! Print to mdout every ntpr steps\n     ntwr= 1,       ! Write a restart file every ntwr steps\n     cut=  10.0,        ! Nonbonded cutoff in Angstroms\n     ntr=1,            ! Turn on restraints\n     restraintmask=":1-314 & !@H=", ! atoms to be restrained\n     restraint_wt=1.0, ! force constant for restraint\n     ntxo=1,           ! Write coordinate file in ASCII format\n     ioutfm=0,         ! Write trajectory file in ASCII format\n     /\n     &qmmm\n     qmmask = \':55\',\n     qmcharge = ,\n     qm_theory = \'XTB\',\n     spin = 2,\n     verbosity = 2,\n     printcharges = 1,\n     qmcut = 10.0,\n     qm_ewald = 1,\n     /\n     &xtb\n     qm_level = \'GFN2-xTB\',\n     tfermi = 300.0,\n     accuracy = 1.d-4\n     maxiter = 250\n     mmhardness = 0.0\n     debug = F\n     /\n    '}
        self._templates.update(sys2_templates)

    @lru_cache(maxsize=8)
    def get_template(self, template_type: str) -> str:
        return self._templates.get(template_type, self._templates['min'])

    def generate_mdin(self, output_path: Path, template_type: str='min') -> str:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(self.get_template(template_type))
        return str(output_path)

    def _generate_files_for_dir(self, output_dir: Path) -> Dict:
        output_dir = Path(output_dir)
        file_specs = [('min.mdin', 'min'), ('min_sys2.mdin', 'min_sys2'), ('sp.mdin', 'sp'), ('vdwm.mdin', 'vdwm'), ('sp_sys2.mdin', 'sp_sys2')]
        try:
            files = {f'{typ}_mdin_file': self.generate_mdin(output_dir / fname, typ) for fname, typ in file_specs}
            return {'output_dir': str(output_dir), **files, 'status': 'success'}
        except Exception as e:
            return {'output_dir': str(output_dir), 'status': 'error', 'error_message': str(e)}

    def generate_batch(self, output_dirs: List[Path]) -> List[Dict]:
        return list(map(self._generate_files_for_dir, output_dirs))

def main():
    generator = MdinGenerator()
if __name__ == '__main__':
    main()