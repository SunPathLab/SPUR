# SPUR Codex UML Diagrams

These diagrams are split by architectural layer to avoid one dense cross-module UML.

## Layered Architecture

```mermaid
flowchart TD
    subgraph ObjectLayer["Object layer"]
        tree["tree.py<br/>TreeNode<br/>tree utilities<br/>constraint_matrix"]
        objects["objects.py<br/>RouteType<br/>GainHistoryTree<br/>Truncal<br/>Branch<br/>Route"]
    end

    subgraph EnumerationLayer["Enumeration layer"]
        enumeration["enumeration.py<br/>TruncalEnumerator<br/>BranchEnumerator<br/>RouteEnumerator"]
    end

    subgraph EncodingLayer["Encoding layer"]
        encoding["encoding.py<br/>RouteEncoder<br/>build_paired_j"]
    end

    subgraph NumericLayer["Scoring / compatibility layer"]
        scoring["scoring.py<br/>TJBScorer<br/>score_tjb<br/>posterior_probabilities<br/>binom_pmf_tensor"]
        compatibility["compatibility.py<br/>PureGainCompatibilityChecker<br/>PureGainFeasibilityCandidate"]
        features["features.py<br/>score_features"]
    end

    subgraph SimulationLayer["Simulation layer"]
        timing["timing.py<br/>sample_timing"]
        simulation["simulation.py<br/>route_jmd<br/>sample_route_jmd"]
        script["scripts/simulate_small_datasets.py<br/>dataset generation"]
        outputs["simulated_small/<br/>metadata.tsv<br/>*.npz"]
    end

    tree --> objects
    tree --> enumeration
    objects --> enumeration
    tree --> encoding
    objects --> encoding
    enumeration --> compatibility
    encoding --> compatibility
    objects --> compatibility
    encoding --> simulation
    objects --> simulation
    tree --> timing
    objects --> timing
    timing --> simulation
    enumeration --> script
    compatibility --> script
    simulation --> script
    scoring --> features
    script --> outputs
```

## Object Model

```mermaid
classDiagram
    direction LR

    class RouteType {
        <<enumeration>>
        ALL
        PURE_GAIN
        MIXED_GAIN_LOSS
    }

    class TreeNode {
        +TreeNode left
        +TreeNode right
        +int id
        +bool is_leaf
    }

    class GainHistoryTree {
        +TreeNode major_tree
        +TreeNode minor_tree
    }

    class Truncal {
        +int major_cn
        +int minor_cn
        +GainHistoryTree tree
        +str tree_hash
        +int n_events
    }

    class Branch {
        +Tuple[int] copy_fate_vector
        +int branch_index
        +bool has_deletion
        +bool is_pure_gain
        +int n_events
    }

    class Route {
        +Truncal truncal
        +Branch branch1
        +Branch branch2
        +bool is_mixed_gain_loss
        +bool is_pure_gain
        +matches_type(RouteType) bool
    }

    TreeNode --> TreeNode : left/right
    GainHistoryTree *-- TreeNode : major_tree
    GainHistoryTree *-- TreeNode : minor_tree
    Truncal *-- GainHistoryTree : tree
    Route *-- Truncal : truncal
    Route *-- Branch : branch1
    Route *-- Branch : branch2
    Route ..> RouteType : filters by
```

## Enumeration Flow

```mermaid
flowchart TD
    input["CN pair input<br/>sample1_cn<br/>sample2_cn<br/>max_mrca_cn<br/>max_branch_events<br/>route_type"] --> routeEnum["RouteEnumerator.enumerate"]

    routeEnum --> truncUpTo["TruncalEnumerator.enumerate_up_to(max_mrca_cn)"]
    truncUpTo --> truncState["enumerate_for_state(major_cn, minor_cn)"]
    truncState --> treeGen["all_full_binary_trees<br/>assign_node_ids<br/>union_tree_hash"]
    treeGen --> truncals["List[Truncal]"]

    truncals --> perTruncal["RouteEnumerator.enumerate_for_truncal"]
    perTruncal --> branch1["BranchEnumerator.enumerate<br/>target = sample1_cn"]
    perTruncal --> branch2["BranchEnumerator.enumerate<br/>target = sample2_cn"]

    branch1 --> allele1["_allele_vectors<br/>major + minor"]
    branch2 --> allele2["_allele_vectors<br/>major + minor"]
    allele1 --> branches1["sample 1 Branch list<br/>filtered by n_events"]
    allele2 --> branches2["sample 2 Branch list<br/>filtered by n_events"]

    branches1 --> product["paired branch cross product"]
    branches2 --> product
    truncals --> product
    product --> route["Route(truncal, branch1, branch2)"]
    route --> filter["Route.matches_type(route_type)"]
    filter --> routes["List[Route]"]
```

## Encoding, Scoring, And Compatibility

```mermaid
flowchart TD
    routes["List[Route]"] --> encoder["RouteEncoder"]

    subgraph Encoding["Encoding"]
        encoder --> copyFate["copy_fate_matrices<br/>branch1, branch2"]
        encoder --> mult["multiplicity_matrices<br/>constraint_matrix.T @ copy_fate.T"]
        mult --> paired["paired_multiplicity_matrices"]
        paired --> jTensor["timing_contribution_tensor<br/>J: routes x timings x joint groups"]
        jTensor --> buildJ["build_paired_j<br/>sample1-major joint index"]
    end

    subgraph Scoring["T/J/B scoring"]
        timingSamples["T timing samples"] --> score["score_tjb"]
        buildJ --> score
        binom["binom_pmf_tensor<br/>B likelihood tensor"] --> score
        score --> logScores["route log scores"]
        logScores --> posterior["posterior_probabilities"]
    end

    subgraph Compatibility["Pure-gain compatibility"]
        q["target JMD q"] --> checker["PureGainCompatibilityChecker.is_compatible"]
        checker --> candidates["_candidates(sample1_cn, sample2_cn)"]
        candidates --> routeEnum["RouteEnumerator<br/>route_type = PURE_GAIN"]
        routeEnum --> candidate["PureGainFeasibilityCandidate"]
        buildJ --> candidate
        constraint["constraint_matrix"] --> candidate
        candidate --> support["support rejection"]
        support --> nnls["NNLS feasibility residual"]
        nnls --> compatible["compatible?"]
    end
```

## Simulation Workflow

```mermaid
flowchart TD
    cnTsv["CN-pair TSV<br/>May_27_real_data_indices.tsv"] --> script["scripts/simulate_small_datasets.py"]
    script --> canonical["canonicalize CN pairs"]
    canonical --> select["select small CN pairs"]

    select --> classes["simulation classes<br/>pure_gain<br/>mixed_gain_loss<br/>conclusive_mixed_gain_loss"]
    classes --> routeEnum["RouteEnumerator.enumerate"]
    routeEnum --> routePool["candidate route pool"]

    routePool --> routePick["sample/select Route"]
    routePick --> sampleRoute["sample_route_jmd"]

    subgraph Timing["Timing sampling"]
        sampleRoute --> sampleTiming["sample_timing(route.truncal, rng)"]
        sampleTiming --> constraints["constraint_matrix"]
        constraints --> feasible["_feasible_timing via linprog"]
        feasible --> randomWalk["null-space random walk<br/>nonnegative timing vector"]
    end

    subgraph JMD["Ground-truth JMD"]
        randomWalk --> routeJmd["route_jmd"]
        routePick --> routeJmd
        routeJmd --> encoder["RouteEncoder.timing_contribution_tensor"]
        encoder --> q["normalized JMD q"]
    end

    q --> compat{"class requires<br/>pure-gain incompatibility?"}
    compat --> checker["PureGainCompatibilityChecker.is_compatible"]
    checker --> keepOrReject["keep/reject example"]
    compat --> keepOrReject

    keepOrReject --> npz["*.npz<br/>jmd<br/>timing<br/>branch copy-fate vectors"]
    keepOrReject --> metadata["metadata.tsv<br/>class label<br/>CN pair<br/>tree hash<br/>branch indices<br/>events<br/>compatibility"]
```
