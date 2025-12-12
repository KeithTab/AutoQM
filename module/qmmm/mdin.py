"""
Optimized MDIN File Generation Module
"""

from pathlib import Path
from typing import List, Dict, Generator
from functools import lru_cache


class MdinGenerator:
    """Optimized MDIN file generator with template caching"""
    
    __slots__ = ('_templates',)
    
    def __init__(self):
        self._templates = {}
        self._init_templates()
    
    def _init_templates(self):
        """Initialize all templates once"""
        self._templates['min'] = """Minmize all the hydrogens
 &cntrl
 ifqnt=1,
 imin=1,           ! Minimize the initial structure
 maxcyc=1,          ! Maximum number of cycles for minimization
 ncyc=1,
 ntb=1,            ! Constant volume
 ntp=0,            ! No pressure scaling
 ntf=1,            ! Complete force evaluation
 ntwx= 1,       ! Write to trajectory file every ntwx steps
 ntpr= 1,       ! Print to mdout every ntpr steps
 ntwr= 1,       ! Write a restart file every ntwr steps
 cut=  10.0,        ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":1-314 & !@H=", ! atoms to be restrained
 restraint_wt=1.0, ! force constant for restraint
 ntxo=1,           ! Write coordinate file in ASCII format
 ioutfm=0,         ! Write trajectory file in ASCII format
 /
 &qmmm
 qmmask = ':111,55',
 qmcharge = 0,
 qm_theory = 'XTB',
 spin = 2,
 verbosity = 2,
 printcharges = 1,
 qmcut = 10.0,
 qm_ewald = 1,
 writepdb = 1,
 /
 &xtb
 qm_level = 'GFN2-xTB',
 tfermi = 300.0,
 accuracy = 1.d-4
 maxiter = 250
 mmhardness = 0.0
 debug = F
 /
"""
        
        self._templates['sp'] = """SP calculation
 &cntrl
 ifqnt=1,
 imin=6,           ! Minimize the initial structure
 nstlim=1,
 dt=0.0,
 ntb=0,            ! Constant volume
 ntp=0,            ! No pressure scaling
 ntf=1,            ! Complete force evaluation
 ntwx= 1,       ! Write to trajectory file every ntwx steps
 ntpr= 1,       ! Print to mdout every ntpr steps
 ntwr= 1,       ! Write a restart file every ntwr steps
 cut=  10.0,        ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":1-314 & !@H=", ! atoms to be restrained
 restraint_wt=1.0, ! force constant for restraint
 ntxo=1,           ! Write coordinate file in ASCII format
 ioutfm=0,         ! Write trajectory file in ASCII format
 /
 &qmmm
 qmmask = ':111,55',
 qmcharge = ,
 qm_theory = 'XTB',
 spin = 2,
 verbosity = 2,
 printcharges = 1,
 qmcut = 10.0,
 qm_ewald = 1,
 writepdb = 1,
 /
 &xtb
 qm_level = 'GFN2-xTB',
 tfermi = 300.0,
 accuracy = 1.d-4
 maxiter = 250
 mmhardness = 0.0
 debug = F
 /
"""
        
        self._templates['vdwm'] = """MMSP calculation
 &cntrl
 imin=6,           ! Minimize the initial structure
 nstlim=1,
 dt=0.0,
 ntb=0,            ! Constant volume
 ntp=0,            ! No pressure scaling
 ntf=1,            ! Complete force evaluation
 ntwx= 1,       ! Write to trajectory file every ntwx steps
 ntpr= 1,       ! Print to mdout every ntpr steps
 ntwr= 1,       ! Write a restart file every ntwr steps
 cut=  10.0,        ! Nonbonded cutoff in Angstroms
 ntr=1,            ! Turn on restraints
 restraintmask=":1-314 & !@H=", ! atoms to be restrained
 restraint_wt=1.0, ! force constant for restraint
 ntxo=1,           ! Write coordinate file in ASCII format
 ioutfm=0,         ! Write trajectory file in ASCII format
 /
"""
        
        sys2_templates = {
            'min_sys2': """Min SYS2 calculation
     &cntrl
     ifqnt=1,
     imin=1,           ! Minimize the initial structure
     maxcyc=1,          ! Maximum number of cycles for minimization
     ncyc=1,
     ntb=1,            ! Constant volume
     ntp=0,            ! No pressure scaling
     ntf=1,            ! Complete force evaluation
     ntwx= 1,       ! Write to trajectory file every ntwx steps
     ntpr= 1,       ! Print to mdout every ntpr steps
     ntwr= 1,       ! Write a restart file every ntwr steps
     cut=  10.0,        ! Nonbonded cutoff in Angstroms
     ntr=1,            ! Turn on restraints
     restraintmask=":1-314 & !@H=", ! atoms to be restrained
     restraint_wt=1.0, ! force constant for restraint
     ntxo=1,           ! Write coordinate file in ASCII format
     ioutfm=0,         ! Write trajectory file in ASCII format
     /
     &qmmm
     qmmask = ':55',
     qmcharge = 0,
     qm_theory = 'XTB',
     spin = 2,
     verbosity = 2,
     printcharges = 1,
     qmcut = 10.0,
     qm_ewald = 1,
     /
     &xtb
     qm_level = 'GFN2-xTB',
     tfermi = 300.0,
     accuracy = 1.d-4
     maxiter = 250
     mmhardness = 0.0
     debug = F
     /
    """,
            'sp_sys2': """SP SYS2 calculation
     &cntrl
     ifqnt=1,
     imin=6,           ! Minimize the initial structure
     nstlim=1,
     dt=0.0,
     ntb=0,            ! Constant volume
     ntp=0,            ! No pressure scaling
     ntf=1,            ! Complete force evaluation
     ntwx= 1,       ! Write to trajectory file every ntwx steps
     ntpr= 1,       ! Print to mdout every ntpr steps
     ntwr= 1,       ! Write a restart file every ntwr steps
     cut=  10.0,        ! Nonbonded cutoff in Angstroms
     ntr=1,            ! Turn on restraints
     restraintmask=":1-314 & !@H=", ! atoms to be restrained
     restraint_wt=1.0, ! force constant for restraint
     ntxo=1,           ! Write coordinate file in ASCII format
     ioutfm=0,         ! Write trajectory file in ASCII format
     /
     &qmmm
     qmmask = ':55',
     qmcharge = ,
     qm_theory = 'XTB',
     spin = 2,
     verbosity = 2,
     printcharges = 1,
     qmcut = 10.0,
     qm_ewald = 1,
     /
     &xtb
     qm_level = 'GFN2-xTB',
     tfermi = 300.0,
     accuracy = 1.d-4
     maxiter = 250
     mmhardness = 0.0
     debug = F
     /
    """
        }
        self._templates.update(sys2_templates)
    
    @lru_cache(maxsize=8)
    def get_template(self, template_type: str) -> str:
        """Get template with caching"""
        return self._templates.get(template_type, self._templates['min'])
    
    def generate_mdin(self, output_path: Path, template_type: str = 'min') -> str:
        """Generate MDIN file from template"""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Use write_text for atomic write
        output_path.write_text(self.get_template(template_type))
        
        return str(output_path)
    
    def _generate_files_for_dir(self, output_dir: Path) -> Dict:
        """Generate all three MDIN files for a directory"""
        output_dir = Path(output_dir)
        
        file_specs = [
            ('min.mdin', 'min'),
            ('min_sys2.mdin', 'min_sys2'),
            ('sp.mdin', 'sp'),
            ('vdwm.mdin', 'vdwm'),
            ('sp_sys2.mdin', 'sp_sys2')
        ]
        
        try:
            # Use dictionary comprehension for file generation
            files = {
                f'{typ}_mdin_file': self.generate_mdin(output_dir / fname, typ)
                for fname, typ in file_specs
            }
            
            return {
                'output_dir': str(output_dir),
                **files,
                'status': 'success'
            }
        except Exception as e:
            return {
                'output_dir': str(output_dir),
                'status': 'error',
                'error_message': str(e)
            }
    
    def generate_batch(self, output_dirs: List[Path]) -> List[Dict]:
        """Generate MDIN files for multiple directories using map"""
        return list(map(self._generate_files_for_dir, output_dirs))


def main():
    generator = MdinGenerator()


if __name__ == '__main__':
    main()
