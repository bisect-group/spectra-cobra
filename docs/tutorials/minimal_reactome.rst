Minimal reactomes
=================

A genome-scale reconstruction is hugely redundant. Most of its reactions are
not needed to grow on any one carbon source, and a large part of the network
can be deleted without the model noticing. The **minimal reactome** is the
smallest set of reactions that still meets a stated requirement.

That question is Burgard, Vaidyaraman and Maranas', from 2001: how few
reactions does *Escherichia coli* need in order to grow? Their answer was
that there is no single answer, because it depends entirely on what is in the
medium — 224 reactions on glucose alone, 229 on acetate alone, and 122 on a
medium engineered to supply everything the cell could transport. A minimal
reactome is minimal **for an environment**, never in the abstract.

This is the other direction from the rest of the package.
:func:`~spectra_cobra.gapfill_for_tasks` adds reactions until a draft can do
something; :func:`~spectra_cobra.minimal_reactome` removes reactions while it
still can.

The problem
-----------

Minimising a *count* is a mixed-integer problem, so this is
:func:`~spectra_cobra.spectra_me` with ``problem_type="minNetMILP"``. There
is no separate solver here: ``minimal_reactome`` builds the problem and hands
it to the same extraction routine everything else in the package uses. Give
every reaction a binary :math:`y_j` tied to its flux,

.. math::

   \min_y \sum_j w_j y_j
   \quad\text{s.t.}\quad
   S v = 0, \quad lb_j\, y_j \le v_j \le ub_j\, y_j

so that :math:`y_j = 0` pins reaction :math:`j` to zero flux and
:math:`\sum_j y_j` is the size of the network. The requirements do not enter
as a core set the way they do in :func:`~spectra_cobra.spectra_me`; they enter
as **bounds**, which is what lets them be a growth rate or a secretion rate
rather than simply "must carry some flux".

There are four kinds of requirement, and you can combine them:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Requirement
     - How it is imposed
   * - ``medium``
     - Every exchange is closed, then the ones named are opened to the
       given uptake rate
   * - ``min_growth`` / ``growth_fraction``
     - A lower bound on the biomass reaction
   * - ``products``
     - A lower bound on the reaction you name, or on the boundary reaction
       the model already has for the metabolite you name
   * - ``tasks``
     - A whole separate condition each — see `Tasks are separate
       conditions`_
   * - several ``medium``\ s
     - A separate condition each, in the same way — see `Several media`_

A first run
-----------

.. code-block:: python

   from cobra.io import load_model
   from spectra_cobra import minimal_reactome

   model = load_model("textbook")
   model.tolerance = 1e-9

   glucose = {"EX_glc__D_e": -10.0, "EX_o2_e": -1000.0, "EX_nh4_e": -1000.0,
              "EX_pi_e": -1000.0, "EX_h2o_e": -1000.0, "EX_h_e": -1000.0,
              "EX_co2_e": -1000.0}

   result = minimal_reactome(model, medium=glucose)
   print(result.summary())

.. code-block:: text

   48 of 95 reactions; growth 0.8739 (needed 0.8739)

``result.model`` is the reduced model, ``result.kept`` and ``result.removed``
the two halves of the reaction list, and ``result.unsatisfied`` is empty —
which is the point of it existing. The result is **verified**: it is rebuilt,
re-solved and re-checked against every requirement, and anything it misses is
reported rather than assumed away.

Nothing is removed from the model you hand in; the reduction happens on a
copy.

The medium decides
------------------

The textbook *E. coli* core model, 95 reactions, with the growth requirement
set to all of what the medium allows and to half of it:

.. list-table::
   :header-rows: 1
   :widths: 30 16 18 18 18

   * - Medium
     - Max growth
     - Reactions at 100%
     - Reactions at 50%
     - Solve
   * - glucose, aerobic
     - 0.8739
     - 48
     - 42
     - < 1 s
   * - acetate, aerobic
     - 0.1733
     - 49
     - 48
     - < 1 s
   * - succinate, aerobic
     - 0.3976
     - 46
     - 45
     - < 1 s
   * - glucose, anaerobic
     - 0.2117
     - 47
     - 43
     - 1 s

Acetate needs *more* reactions than glucose despite supporting a fifth of the
growth. Comparing the two sets at 100% says why: the acetate network has to
carry the glyoxylate shunt (``ICL``, ``MALS``), gluconeogenesis (``FBP``,
``PPCK``) and acetate activation (``ACKr``, ``PTAr``), none of which glucose
needs, while glucose gets by with glycolysis and the pentose phosphate
pathway (``GLCpts``, ``PFK``, ``PYK``, ``PDH``, ``G6PDH2r``, ``GND``,
``PGL``, ``PPC``).

.. code-block:: python

   set(acetate.kept) - set(glucose.kept)

.. code-block:: text

   {'ACKr', 'ACt2r', 'EX_ac_e', 'FBP', 'ICL', 'MALS',
    'ME1', 'NADTRHD', 'PPCK', 'PTAr'}

This is Burgard's finding reproduced in miniature — the size of a minimal
reactome tracks the difficulty of the medium, not the growth rate it
supports.

Several media
-------------

Pass more than one medium and the result is the smallest network that works
in **all** of them — not the smaller of the two answers. Each medium becomes
its own condition with its own flux vector, exactly as a task does, and
``growth_fraction`` is measured separately in each.

.. code-block:: python

   result = minimal_reactome(
       model,
       medium={"glucose": glucose, "acetate": acetate, "anaerobic": anaerobic},
       growth_fraction=1.0,
   )
   print(result.conditions)
   print({k: round(v, 4) for k, v in result.growth.items()})

.. code-block:: text

   ('glucose', 'acetate', 'anaerobic')
   {'glucose': 0.8739, 'acetate': 0.1733, 'anaerobic': 0.2117}

On the textbook model those three alone need 48, 49 and 47 reactions, whose
union is 66; solved jointly the answer is **65**. The same thing on iJO1366,
with glucose, acetate and glycerol at full growth:

.. list-table::
   :header-rows: 1
   :widths: 60 40

   * -
     - Reactions
   * - glucose alone
     - 437
   * - acetate alone
     - 431
   * - glycerol alone
     - 433
   * - **union of the three**
     - **458**
   * - **solved jointly** (24 s)
     - **448**

Ten reactions, and 415 of them are common to all three answers anyway. This
is NetworkReducer's original demonstration — a reduced network that keeps aerobic *and*
anaerobic growth — and the reason to do it in one solve is the same reason
as for tasks: separate answers each pick whatever route is cheapest for
them, and nothing pushes them onto shared routes.

A plain list works too and is labelled for you (``medium_0``, ``medium_1``,
…); a mapping lets you name them. Because ``growth``, ``required_growth``
and ``production`` are then keyed by label, they are keyed by label in the
single-medium case too, under the label ``"growth"``:

.. code-block:: python

   result.growth["growth"]          # one medium
   result.growth["acetate"]         # several

How much growth to keep
-----------------------

Two mutually exclusive ways to say it:

``growth_fraction``
   A fraction of what the model achieves **in this medium** before anything is
   removed. The default is ``1.0`` — lose no growth at all, which is MinReact's
   own default. The achievable rate is measured first, which is why asking for
   a fraction of a model that cannot grow in the medium is an error rather
   than a silent zero.

``min_growth``
   An absolute lower bound on the biomass reaction, in the model's own units.

.. code-block:: python

   minimal_reactome(model, medium=glucose, growth_fraction=0.5)
   minimal_reactome(model, medium=glucose, min_growth=0.4)

Pass ``growth_fraction=0`` to not require growth at all, which is what you
want when the requirement is a task list and growth is beside the point.

Requiring a product
-------------------

``products`` maps a reaction or a metabolite to the flux it must carry.
Naming a **reaction** puts a lower bound on it. Naming a **metabolite**
resolves to the boundary reaction the model already has for it, so
``products={"ac_e": 1.0}`` and ``products={"EX_ac_e": 1.0}`` are the same
request.

.. important::

   **Nothing is added to the network.** A metabolite the model does not
   exchange is refused:

   .. code-block:: text

      'y_c' has no exchange reaction, so its production cannot be required
      without adding one -- and adding one would change the network being
      minimised.

   Quietly adding a demand reaction would be the obvious convenience and it
   is the wrong answer. The result would no longer be a subnetwork of *your*
   model but of a model you do not have, and it would be feasible only
   because of the drain that was added to measure it — take the drain away
   and the same reaction set may produce nothing at all. If an internal
   metabolite really is the thing you want, add the reaction to your own
   model first, where it is visible and counted like everything else.

   Two boundary reactions for the same metabolite is also refused, since
   which one the production is measured through would otherwise be
   arbitrary, and so is an exchange written as a *source* (``--> met``),
   where flux is uptake and a lower bound would require the opposite of
   production.

.. code-block:: python

   result = minimal_reactome(model, medium=glucose, growth_fraction=0.5,
                             products={"EX_succ_e": 5.0})

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 30

   * - Succinate demanded
     - Reactions
     - Growth reached
     - Succinate made
   * - none
     - 42
     - 0.7040
     - —
   * - ≥ 1.0
     - 45
     - 0.5744
     - 5.635
   * - ≥ 5.0
     - 45
     - 0.5744
     - 5.635
   * - ≥ 8.0
     - 48
     - 0.4400
     - 8.000

.. important::

   Growth and production are **one condition**, solved as one flux vector.
   A network that can grow *or* secrete, but not both at once, has not met
   the requirement. Demanding 10 units of succinate alongside half the
   maximum growth is impossible on this medium, and it is refused before the
   MILP is built:

   .. code-block:: text

      The model cannot meet the growth condition even with every reaction
      available, so no subnetwork of it can either. Required together in one
      flux distribution: Biomass_Ecoli_core at 0.436961 (reachable);
      EX_succ_e at 10 (reachable). Growth and production share a medium, so
      each can be reachable on its own and the pair still not be. Lower one
      of them, or open the medium.

   Each is reachable alone. Together they are not. Checking the full network
   first is what turns an infeasible MILP — which says nothing about *which*
   requirement was impossible — into an answer.

Tasks are separate conditions
-----------------------------

A metabolic task is not another bound on the same flux vector.
:func:`~spectra_cobra.tasks.task_constraints` closes the model's own
exchanges and opens its own, so every task is a **different environment**, and
a network that performs two tasks does so with two different flux
distributions. The reactions are shared; the fluxes are not.

So the problem becomes one flux vector per condition over a single set of
binaries:

.. math::

   \min_y \sum_j w_j y_j
   \quad\text{s.t.}\quad
   S^{(k)} v^{(k)} = 0, \quad
   lb^{(k)}_j\, y_j \le v^{(k)}_j \le ub^{(k)}_j\, y_j
   \quad \forall k

This is the "protected functions" idea of NetworkReducer and MinNW: a reduced
network must preserve not one phenotype but a list of them. Internally the
network is replicated, one copy per condition, with no metabolite shared
between the copies — so the model handed to the solver is *k* times the size
of yours, and a long task list is what makes this expensive.

The replicated model then goes to :func:`~spectra_cobra.spectra_me`, with
each reaction's copies declared as one ``indicator_groups`` entry so that a
single binary governs all of them. The whole multi-condition machinery is
that one argument; everything else is the ordinary extraction.

.. code-block:: python

   from spectra_cobra import MetabolicTask, TaskEquation, minimal_reactome

   hydrolysis = {"atp_c": -1.0, "h2o_c": -1.0,
                 "adp_c": 1.0, "pi_c": 1.0, "h_c": 1.0}
   aerobic = ["glc__D_e", "o2_e", "nh4_e", "pi_e", "h2o_e", "h_e"]
   outputs = {m: (0.0, 1000.0) for m in
              ["co2_e", "h2o_e", "h_e", "ac_e", "etoh_e",
               "for_e", "lac__D_e", "succ_e"]}

   tasks = [
       MetabolicTask(
           id="atp_aerobic_glucose",
           inputs={m: (0.0, 10.0) for m in aerobic},
           outputs=outputs,
           equations=[TaskEquation(hydrolysis, (1.0, 1000.0))],
       ),
       MetabolicTask(
           id="atp_anaerobic_glucose",
           inputs={m: (0.0, 10.0) for m in aerobic if m != "o2_e"},
           outputs=outputs,
           equations=[TaskEquation(hydrolysis, (1.0, 1000.0))],
       ),
       MetabolicTask(
           id="akg_from_acetate",
           inputs={m: (0.0, 10.0) for m in
                   ["ac_e", "o2_e", "nh4_e", "pi_e", "h2o_e", "h_e"]},
           outputs=outputs,
           equations=[TaskEquation({"akg_c": -1.0}, (1.0, 1000.0))],
       ),
   ]

   result = minimal_reactome(model, medium=glucose, growth_fraction=0.5,
                             tasks=tasks)
   print(result.conditions, len(result.kept), result.unsatisfied)

.. code-block:: text

   ('growth', 'atp_aerobic_glucose', 'atp_anaerobic_glucose',
    'akg_from_acetate') 51 ()

Why solve them jointly
~~~~~~~~~~~~~~~~~~~~~~

Because the union of separate answers is not the smallest network that
satisfies them all. On the textbook model, with the four conditions above:

.. list-table::
   :header-rows: 1
   :widths: 60 20 20

   * - Condition
     - Alone
     - Jointly
   * - growth on glucose at 50%
     - 42
     -
   * - ``atp_aerobic_glucose``
     - 13
     -
   * - ``atp_anaerobic_glucose``
     - 13
     -
   * - ``akg_from_acetate``
     - 22
     -
   * - **union of the four**
     - **58**
     - **51**

Seven reactions, on a 95-reaction model. Each condition, solved on its own,
picks whatever route is cheapest *for it*; solved together they are pushed
onto shared routes that no one of them would have chosen. That is the whole
reason the conditions share one set of binaries, and it is the difference
between this and running :func:`~spectra_cobra.minimal_reactome` once per
condition and taking the union.

.. note::

   A task the **full** model cannot perform is refused up front, by name.
   Removing reactions cannot make a task work, so such a task makes the MILP
   infeasible, and an infeasible MILP does not say which task was at fault.
   If you want those tasks, gap-fill towards them first with
   :func:`~spectra_cobra.gapfill_for_tasks`.

Keeping reactions, and making them work
---------------------------------------

``keep_reactions`` pins reactions into the answer, and it does more than
retain them. Each one **is** the core set handed to
:func:`~spectra_cobra.spectra_me`: the two LPs of its direction phase settle
which way the reaction runs, and the bound rewrite then forces at least
``tol`` of flux through it — the same mechanism, and the same code, as a core
reaction in any other extraction.

.. code-block:: python

   minimal_reactome(model, medium=glucose,
                    keep_reactions=["ATPM", "CYTBD"], tol=1e-4)

That distinction matters, because retaining a reaction is not the same as the
reaction working. A reaction kept but sitting at zero flux is dead weight in
the result. On the toy cell of the test suite, keeping the acetate transporter
``T_AC`` on a glucose medium pulls ``P_AC`` in as well — ``P_AC`` is the only
thing that can make the acetate ``T_AC`` consumes, and nothing else on that
medium takes acetate up:

.. code-block:: text

   keep_reactions=["T_AC"]            -> EX_glc_e, T_GLC, BIO, T_AC, P_AC
   keep_reactions=["T_GLC_LONG_1"]    -> EX_glc_e, T_GLC_LONG_1, T_GLC_LONG_2, BIO

Merely retaining them would have left both ``T_AC`` and ``T_GLC_LONG_1`` at
zero and neither partner in the model.

.. important::

   The flux requirement applies to the **growth conditions only**. A task
   supplies its own medium and closes every exchange the model has, so forcing
   your reactions inside someone else's environment would be surprising — and
   for any boundary reaction, infeasible by construction: keeping
   ``EX_glc_e`` alongside any task at all would fail. Inside a task a kept
   reaction is simply present.

A reaction that cannot reach ``tol`` in some growth condition is retained
anyway, and said so rather than silently forced or silently dropped:

.. code-block:: python

   result.kept_but_blocked

.. code-block:: text

   {'growth': ('DEAD',)}

with a warning naming them. Note that ``keep_reactions`` is exempt from
preprocessing in both directions: such a reaction is never deleted as blocked,
and no kept reaction is tested for essentiality, since it is kept either way.

Weights
-------

``weights`` is what each reaction *costs* to keep — not how much flux it gets.
The default is 1 each, so the objective is a plain count; raising one makes
the solver prefer a longer route around it, which is how evidence enters if
you have any.

.. code-block:: python

   minimal_reactome(model, medium=glucose,
                    weights={"PGI": 5.0})   # prefer any route but this one

Weights must be non-negative: a negative cost would pay the solver to keep
reactions. A weight of 0 leaves a reaction free but still removable — use
``keep_reactions`` to force one in. A weight on a reaction that has no binary,
because it was kept or proved essential, does nothing and is reported.

Alternative minimal reactomes
-----------------------------

A redundant network has many equally small answers, so treat one run as *a*
minimal reactome rather than *the* minimal reactome. ``n_solutions`` forbids
each answer from recurring before solving again:

.. code-block:: python

   results = minimal_reactome(model, medium=glucose, growth_fraction=0.5,
                              n_solutions=5)
   print([len(r.kept) for r in results])
   print(len(set.intersection(*[set(r.kept) for r in results])))

.. code-block:: text

   [42, 42, 43, 43, 44]
   34

Only the first is guaranteed minimal; the exclusion is a no-good cut on the
binary assignment, not an enumeration of the optima, so later answers may be
larger. The 34 reactions common to all five are the more interesting number:
those are the ones no alternative route exists for. If fewer alternatives
exist than you asked for, you get what there is, with a warning.

At genome scale
---------------

*E. coli* iJO1366 — 2583 reactions, 1805 metabolites — at **100%** of the
growth each medium allows, so nothing is given up at all:

.. list-table::
   :header-rows: 1
   :widths: 20 16 14 12 18 14

   * - Medium
     - Max growth
     - Reactions
     - Of 2583
     - Of those, essential
     - Solve
   * - glucose
     - 0.98237
     - 437
     - 17%
     - 407
     - 14 s
   * - acetate
     - 0.24720
     - 431
     - 17%
     - 400
     - 15 s
   * - glycerol
     - 0.56279
     - 433
     - 17%
     - 400
     - 14 s
   * - succinate
     - 0.49248
     - 434
     - 17%
     - 401
     - 15 s

Every one is verified: the reduced model rebuilt from scratch reaches the full
wild-type rate to six figures. Six sevenths of iJO1366 is doing nothing on a
glucose minimal medium — and of the seventh that is left, more than nine
tenths has no alternative route at all.

These are in the same neighbourhood as the published figures — MinReact
reports 430 reactions for an *E. coli* K-12 model against MinNW's 447 — but
that is not a like-for-like comparison: different reconstruction, different
medium, different growth cutoff. Burgard's 224 was a 2001-era flux balance
model an order of magnitude smaller than iJO1366, so the two numbers are not
comparable either.

.. warning::

   **Full growth is the easy case.** Demanding 100% of the optimum leaves the
   solver very little freedom and the MILP closes in about fifteen seconds.
   Relaxing to 50% opens up an enormous space of equally feasible networks,
   and the same problem does not close at all — it was still unproven at a
   1200 s limit.

   What comes back then is the incumbent, with a warning that it is feasible
   but may not be minimal. It is usually a perfectly good answer to "a small
   network that works". At an equal 1200 s budget on glucose at 50%:

   .. list-table::
      :header-rows: 1
      :widths: 40 20 20 20

      * -
        - Reactions
        - Essential
        - Verified
      * - ``preprocess=True``
        - 425
        - 295
        - yes
      * - ``preprocess=False``
        - 440
        - —
        - yes

   But not always. Before preprocessing existed, an 1800 s run of this same
   problem returned a network that **failed** verification, with
   ``result.unsatisfied == ('growth',)``, while a 120 s run of it returned
   440 reactions that verified. The longer search got closer to the true
   optimum and in doing so leaned on fluxes too small to survive the rebuild.

   This is why ``check`` defaults to on and why ``unsatisfied`` exists. Test
   it, rather than assuming a returned model is a valid one:

   .. code-block:: python

      result = minimal_reactome(model, medium=glucose, growth_fraction=0.5)
      if result.unsatisfied:
          raise RuntimeError(f"does not meet {result.unsatisfied}")

   If a solve stops at its limit and fails, give it more time so it can
   prove optimality, or lower ``inclusion_cutoff`` so fewer trace-carrying
   reactions are discarded.

Deciding what can be decided first
----------------------------------

``preprocess`` is on by default, and it takes two exact decisions before the
MILP starts.

**Blocked.** A reaction that can carry no flux in *any* condition is in no
feasible solution of any of them, so it is deleted outright — which shrinks
every one of the *k* LP copies, not just the binary count. This is
:func:`~spectra_cobra.spectra_cc` run inside each condition and
**intersected**: a reaction useless under the growth medium may be the only
route a task has.

**Essential.** A reaction some condition is *infeasible without* is in every
feasible solution, so it is forced in and given no binary. **Unioned** over
conditions, for the mirror-image reason. The test is feasibility, never a
growth comparison: because the growth condition already carries the biomass
lower bound, "infeasible without this reaction" *is* "cannot reach the
required growth rate without it" — and that same test transfers unchanged to
a task, which has no growth at all.

On iJO1366 on glucose minimal medium, at 100% of wild-type growth:

.. list-table::
   :header-rows: 1
   :widths: 60 20 20

   * -
     - Reactions
     - Share
   * - blocked in every condition, deleted
     - 2057
     - 80%
   * - essential somewhere, no binary
     - 407
     - 16%
   * - **binaries left to decide**
     - **119**
     - **5%**

Ninety-five per cent of the integer problem is settled before the solver
sees it. Both passes are exact, and the objective confirms it: **28 either
way**, with ``preprocess=False`` reporting 435 reactions and
``preprocess=True`` 437 — the same optimum, read off two solutions whose
trace fluxes happen to sit in different places.

.. note::

   It is not free, and on an easy problem it does not pay. The same glucose
   run takes 7 s with ``preprocess=False`` and 14 s with it on — the
   consistency checks and knockouts cost more than the MILP they shorten.
   Turn it off for small models and for problems that already close
   quickly; leave it on when the MILP is the bottleneck.

   Where it does pay is the problem that does not close. On glucose at 50%
   growth, at an equal 1200 s budget, it returns 425 reactions against 440
   without it — fifteen fewer, for the same wait.

   Note also that the reported network can differ by a reaction or two
   between the two settings — 435 against 437 at full growth — even though
   the optimum is identical. That is the trace-flux rule below, not a
   disagreement about what is minimal.

``result.essential`` carries the reactions proved to be in every possible
answer. That is arguably the more interesting output than the size: they are
the part of the network with no alternative route at all.

Two things to watch
-------------------

**The binaries are not the whole answer.** A binary the solver leaves at its
integrality tolerance rather than at exactly zero still multiplies a flux
bound of up to a thousand, so a reaction the solution genuinely depends on can
be reported as switched off. Reading membership off the binaries alone
produced, on iJO1366 on acetate, a 412-reaction network that could not grow at
all — one of the leaked reactions was the iron exchange. Any reaction carrying
more than ``inclusion_cutoff`` (1e-11 by default) is therefore kept whatever
its binary says, which brought that case to 434 reactions and full growth.
The count can consequently be a little above the MILP's own objective value.
That trade is deliberate, and it is the same one
:func:`~spectra_cobra.spectra_me` makes.

**Set** ``model.tolerance`` **before you start.** The default of 1e-7 is loose
enough on a genome-scale model for a requirement to be met only to within it.
Use 1e-9, the lowest most solvers accept.

.. code-block:: python

   model.tolerance = 1e-9

Next
----

* :doc:`../functions/reactome` — the API
* :doc:`formulations` — what ``minNetMILP`` is doing, and the alternatives
* :doc:`gapfilling_tasks` — the opposite operation, with the same tasks
* :doc:`minimal_microbiome` — the same question asked of organisms

References
----------

* Burgard, A. P., Vaidyaraman, S., and Maranas, C. D. (2001). Minimal
  reaction sets for *Escherichia coli* metabolism under different growth
  requirements and uptake environments. *Biotechnology Progress*, 17(5),
  791–797. https://doi.org/10.1021/bp0100880
* Sambamoorthy, G., and Raman, K. (2020). MinReact: a systematic approach for
  identifying minimal metabolic networks. *Bioinformatics*, 36(15),
  4309–4315. https://doi.org/10.1093/bioinformatics/btaa497
* Erdrich, P., Steuer, R., and Klamt, S. (2015). An algorithm for the
  reduction of genome-scale metabolic network models to meaningful core
  models. *BMC Systems Biology*, 9, 48.
  https://doi.org/10.1186/s12918-015-0191-x
* Röhl, A., and Bockmayr, A. (2017). A mixed-integer linear programming
  approach to the reduction of genome-scale metabolic networks. *BMC
  Bioinformatics*, 18, 2. https://doi.org/10.1186/s12859-016-1412-z
* Jonnalagadda, S., and Srinivasan, R. (2014). An efficient graph theory
  based method to identify every minimal reaction set in a metabolic network.
  *BMC Systems Biology*, 8, 28. https://doi.org/10.1186/1752-0509-8-28
