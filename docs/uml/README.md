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
    input["CN pair input<br/>sample1_cn<br/>sample2_cn<br/>max_truncal_cn<br/>max_branch_events<br/>route_type"] --> routeEnum["RouteEnumerator.enumerate"]

    routeEnum --> truncUpTo["TruncalEnumerator.enumerate_up_to(max_truncal_cn)"]
    truncUpTo --> truncState["enumerate_for_state(major_cn, minor_cn)"]
    truncState --> treeGen["all_full_binary_trees<br/>assign_node_ids<br/>union_tree_shape"]
    treeGen --> truncals["List[Truncal]"]

    truncals --> perTruncal["RouteEnumerator.enumerate_for_truncal"]
    perTruncal --> branch1["BranchEnumerator.enumerate<br/>target = sample1_cn"]
    perTruncal --> branch2["BranchEnumerator.enumerate<br/>target = sample2_cn"]

    branch1 --> allele1["_allele_vectors<br/>major + minor"]
    branch2 --> allele2["_allele_vectors<br/>major + minor"]
    allele1 --> branchCap1["enumeration cap<br/>branch1.n_events <= max_branch_events"]
    allele2 --> branchCap2["enumeration cap<br/>branch2.n_events <= max_branch_events"]
    branchCap1 --> branches1["sample 1 Branch list"]
    branchCap2 --> branches2["sample 2 Branch list"]

    branches1 --> product["paired branch cross product"]
    branches2 --> product
    truncals --> product
    product --> route["Route(truncal, branch1, branch2)"]
    route --> typeFilter["Route.matches_type(route_type)"]
    typeFilter --> routes["List[Route]"]
    routes --> postFilter{"optional post-enumeration<br/>route_filter"}
    postFilter --> noFilter["none"]
    postFilter --> total["total_events <= threshold"]
    postFilter --> maxBranch["max(branch1, branch2) events <= threshold"]
    postFilter --> branchSum["branch1 + branch2 events <= threshold"]
    postFilter --> singleSample["max(truncal + branch events) <= threshold"]
    noFilter --> finalRoutes["filtered Route set"]
    total --> finalRoutes
    maxBranch --> finalRoutes
    branchSum --> finalRoutes
    singleSample --> finalRoutes
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
    input["CN-pair input<br/>--cn-pair total:minor_total:minor"] --> script["scripts/simulate_datasets.py"]
    script --> parse["parse_between_sample_cn"]
    parse --> truncalCap["resolve max_truncal_cn<br/>min(max total CN + 2, 6) unless provided"]

    truncalCap --> routeEnum["RouteEnumerator.enumerate"]
    routeEnum --> routeType["split by route type<br/>pure_gain or mixed_gain_loss"]
    routeType --> routeFilter{"optional route_filter"}
    routeFilter --> routePool["filtered route pool"]

    routePool --> perRoute["for each Route"]

    subgraph Timing["Timing sampling"]
        perRoute --> sampleTiming["TimingSampler.sample_timing_vectors"]
        sampleTiming --> constraints["route.truncal.constraint_matrix"]
        constraints --> feasible["_feasible_timing via NNLS"]
        feasible --> hitAndRun["hit_and_run<br/>nonnegative timing vectors T"]
    end

    subgraph JMD["Ground-truth JMD"]
        perRoute --> h["route.h(dtype)"]
        hitAndRun --> th["Q = T @ H"]
        h --> th
        th --> q["normalized JMD rows"]
    end

    q --> compat{"route_type is<br/>conclusive_mixed_gain_loss?"}
    compat --> checker["PureGainCompatibilityChecker.is_compatible"]
    checker --> keepOrReject["drop pure-gain-compatible Q rows"]
    compat --> keepOrReject
    keepOrReject --> densityFilter["density threshold filter"]
    densityFilter --> sampleQ["sample Q rows"]
    sampleQ --> vaf["JointVAFSampler.sample_joint_vafs"]
    vaf --> reads["generate_total_reads<br/>generate_alt_reads"]

    reads --> csv["read-count CSV files"]
    sampleQ --> truth["ground-truth density"]
    vaf --> truthVaf["ground-truth expected VAFs"]
    script --> metadata["metadata.json"]
```
