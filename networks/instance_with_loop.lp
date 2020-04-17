reaction("importS1").
product("S1", "importS1").

reaction("r1").
reactant("S1", "r1").
reactant("I2", "r1").
product("I1", "r1").

reaction("r2").
reactant("I1", "r2").
product("I2", "r2").
product("T", "r2").