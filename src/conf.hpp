/*
  tkrsadi
  Turbocharger parameters calculation.

  File: conf.hpp

  Copyright (C) 2021 Artem Petrov <pa23666@yandex.ru>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CONF_HPP
#define CONF_HPP

class Conf {

public:

    bool readConfigFile();

    double val_Ne()       const { return m_Ne;       }
    double val_ge()       const { return m_ge;       }
    double val_alpha()    const { return m_alpha;    }
    double val_phip()     const { return m_phip;     }
    double val_pik()      const { return m_pik;      }
    double val_itkr()     const { return m_itkr;     }
    double val_p0()       const { return m_p0;       }
    double val_t0()       const { return m_t0;       }
    double val_sigma_in() const { return m_sigma_in; }
    double val_ca()       const { return m_ca;       }
    double val_etakn()    const { return m_etakn;    }
    double val_phi()      const { return m_phi;      }
    double val_zeta_in()  const { return m_zeta_in;  }
    double val_D0D1()     const { return m_D0D1;     }
    double val_D1D2()     const { return m_D1D2;     }
    double val_ia()       const { return m_ia;       }
    double val_tau1()     const { return m_tau1;     }
    double val_zeta1()    const { return m_zeta1;    }
    double val_zeta2()    const { return m_zeta2;    }
    double val_kc2r()     const { return m_kc2r;     }
    double val_alphap()   const { return m_alphap;   }
    double val_zk()       const { return m_zk;       }
    double val_tau2()     const { return m_tau2;     }
    double val_kb3()      const { return m_kb3;      }
    double val_kD3()      const { return m_kD3;      }
    double val_zetad()    const { return m_zetad;    }
    double val_kD4()      const { return m_kD4;      }
    double val_kb4()      const { return m_kb4;      }
    double val_kalpha4()  const { return m_kalpha4;  }
    double val_zd()       const { return m_zd;       }
    double val_n4()       const { return m_n4;       }
    double val_tau3()     const { return m_tau3;     }
    double val_tau4()     const { return m_tau4;     }
    double val_nul()      const { return m_nul;      }
    double val_kck()      const { return m_kck;      }

    double val_pr_z()     const { return m_pr_z;     }
    double val_tr_z()     const { return m_tr_z;     }
    double val_sigmay()   const { return m_sigmay;   }
    double val_kD1t()     const { return m_kD1t;     }
    double val_etat()     const { return m_etat;     }
    double val_kD2t()     const { return m_kD2t;     }
    double val_kDbt()     const { return m_kDbt;     }
    double val_rhot()     const { return m_rhot;     }
    double val_phit()     const { return m_phit;     }
    double val_alpha1t()  const { return m_alpha1t;  }
    double val_psi()      const { return m_psi;      }
    double val_delta()    const { return m_delta;    }
    double val_kNrd()     const { return m_kNrd;     }
    double val_etamt()    const { return m_etamt;    }

private:

    bool createBlank() const;

    double m_Ne       = 0;
    double m_ge       = 0;
    double m_alpha    = 0;
    double m_phip     = 0;
    double m_pik      = 0;
    double m_itkr     = 0;
    double m_p0       = 0;
    double m_t0       = 0;
    double m_sigma_in = 0;
    double m_ca       = 0;
    double m_etakn    = 0;
    double m_phi      = 0;
    double m_zeta_in  = 0;
    double m_D0D1     = 0;
    double m_D1D2     = 0;
    double m_ia       = 0;
    double m_tau1     = 0;
    double m_zeta1    = 0;
    double m_zeta2    = 0;
    double m_kc2r     = 0;
    double m_alphap   = 0;
    double m_zk       = 0;
    double m_tau2     = 0;
    double m_kb3      = 0;
    double m_kD3      = 0;
    double m_zetad    = 0;
    double m_kD4      = 0;
    double m_kb4      = 0;
    double m_kalpha4  = 0;
    double m_zd       = 0;
    double m_n4       = 0;
    double m_tau3     = 0;
    double m_tau4     = 0;
    double m_nul      = 0;
    double m_kck      = 0;

    double m_pr_z     = 0;
    double m_tr_z     = 0;
    double m_sigmay   = 0;
    double m_kD1t     = 0;
    double m_etat     = 0;
    double m_kD2t     = 0;
    double m_kDbt     = 0;
    double m_rhot     = 0;
    double m_phit     = 0;
    double m_alpha1t  = 0;
    double m_psi      = 0;
    double m_delta    = 0;
    double m_kNrd     = 0;
    double m_etamt    = 0;

};

#endif // CONF_HPP
