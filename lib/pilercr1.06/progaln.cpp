#include "multaln.h"

#define TRACE		0
#define VALIDATE	0
#define TRACE_LENGTH_DELTA	0

static bool g_bVerbose = false;

static void LogLeafNames(const Tree &tree, unsigned uNodeIndex)
	{
	const unsigned uNodeCount = tree.GetNodeCount();
	unsigned *Leaves = new unsigned[uNodeCount];
	unsigned uLeafCount;
	GetLeaves(tree, uNodeIndex, Leaves, &uLeafCount);
	for (unsigned i = 0; i < uLeafCount; ++i)
		{
		if (i > 0)
			Log(",");
		Log("%s", tree.GetLeafName(Leaves[i]));
		}
	delete[] Leaves;
	}

void ProgressiveAlign(const SeqVect &Seqs, const Tree &GuideTree, MSA &Aln)
	{
	assert(Seqs.GetSeqCount() > 0);
	assert(GuideTree.IsRooted());

#if	TRACE
	Log("GuideTree:\n");
	GuideTree.LogMe();
#endif

	const unsigned uSeqCount = Seqs.Length();
	const unsigned uNodeCount = 2*uSeqCount - 1;
	const unsigned uIterCount = uSeqCount - 1;

	ProgNode *ProgNodes = new ProgNode[uNodeCount];

	unsigned uJoin = 0;
	unsigned uTreeNodeIndex = GuideTree.FirstDepthFirstNode();
	do
		{
		if (GuideTree.IsLeaf(uTreeNodeIndex))
			{
			if (uTreeNodeIndex >= uNodeCount)
				Quit("TreeNodeIndex=%u NodeCount=%u\n", uTreeNodeIndex, uNodeCount);
			ProgNode &Node = ProgNodes[uTreeNodeIndex];
			unsigned uId = GuideTree.GetLeafId(uTreeNodeIndex);
			if (uId >= uSeqCount)
				Quit("Seq index out of range");
			const Seq &s = *(Seqs[uId]);
			Node.m_MSA.FromSeq(s);
			Node.m_MSA.SetSeqId(0, uId);
			Node.m_uLength = Node.m_MSA.GetColCount();
		// TODO: Term gaps settable
			Node.m_Prof = ProfileFromMSA(Node.m_MSA);
			Node.m_EstringL = 0;
			Node.m_EstringR = 0;
#if	TRACE
			Log("Leaf id=%u\n", uId);
			Log("MSA=\n");
			Node.m_MSA.LogMe();
			Log("Profile (from MSA)=\n");
			ListProfile(Node.m_Prof, Node.m_uLength, &Node.m_MSA);
#endif
			}
		else
			{
			++uJoin;

			const unsigned uMergeNodeIndex = uTreeNodeIndex;
			ProgNode &Parent = ProgNodes[uMergeNodeIndex];

			const unsigned uLeft = GuideTree.GetLeft(uTreeNodeIndex);
			const unsigned uRight = GuideTree.GetRight(uTreeNodeIndex);

			if (g_bVerbose)
				{
				Log("Align: (");
				LogLeafNames(GuideTree, uLeft);
				Log(") (");
				LogLeafNames(GuideTree, uRight);
				Log(")\n");
				}

			ProgNode &Node1 = ProgNodes[uLeft];
			ProgNode &Node2 = ProgNodes[uRight];

#if	TRACE
			Log("AlignTwoMSAs:\n");
#endif
			AlignProfiles(Node1.m_Prof, Node1.m_uLength, Node2.m_Prof,
			  Node2.m_uLength, Parent.m_Path, &Parent.m_Prof, &Parent.m_uLength);

#if	TRACE_LENGTH_DELTA
			{
			unsigned L = Node1.m_uLength;
			unsigned R = Node2.m_uLength;
			unsigned P = Parent.m_Path.GetEdgeCount();
			unsigned Max = L > R ? L : R;
			unsigned d = P - Max;
			Log("LD%u;%u;%u;%u\n", L, R, P, d);
			}
#endif
			PathToEstrings(Parent.m_Path, &Parent.m_EstringL, &Parent.m_EstringR);

#if	VALIDATE
			{
#if	TRACE
			Log("AlignTwoMSAs:\n");
#endif
			PWPath TmpPath;
			AlignTwoMSAs(Node1.m_MSA, Node2.m_MSA, Parent.m_MSA, TmpPath);
			ProfPos *P1 = ProfileFromMSA(Node1.m_MSA, true);
			ProfPos *P2 = ProfileFromMSA(Node2.m_MSA, true);
			unsigned uLength = Parent.m_MSA.GetColCount();
			ProfPos *TmpProf = ProfileFromMSA(Parent.m_MSA, true);

#if	TRACE
			Log("Node1 MSA=\n");
			Node1.m_MSA.LogMe();

			Log("Node1 prof=\n");
			ListProfile(Node1.m_Prof, Node1.m_MSA.GetColCount(), &Node1.m_MSA);
			Log("Node1 prof (from MSA)=\n");
			ListProfile(P1, Node1.m_MSA.GetColCount(), &Node1.m_MSA);

			AssertProfsEq(Node1.m_Prof, Node1.m_uLength, P1, Node1.m_MSA.GetColCount());

			Log("Node2 prof=\n");
			ListProfile(Node2.m_Prof, Node2.m_MSA.GetColCount(), &Node2.m_MSA);

			Log("Node2 MSA=\n");
			Node2.m_MSA.LogMe();

			Log("Node2 prof (from MSA)=\n");
			ListProfile(P2, Node2.m_MSA.GetColCount(), &Node2.m_MSA);

			AssertProfsEq(Node2.m_Prof, Node2.m_uLength, P2, Node2.m_MSA.GetColCount());

			TmpPath.AssertEqual(Parent.m_Path);

			Log("Parent MSA=\n");
			Parent.m_MSA.LogMe();

			Log("Parent prof=\n");
			ListProfile(Parent.m_Prof, Parent.m_uLength, &Parent.m_MSA);

			Log("Parent prof (from MSA)=\n");
			ListProfile(TmpProf, Parent.m_MSA.GetColCount(), &Parent.m_MSA);

#endif	// TRACE
			AssertProfsEq(Parent.m_Prof, Parent.m_uLength,
			  TmpProf, Parent.m_MSA.GetColCount());
			delete[] P1;
			delete[] P2;
			delete[] TmpProf;
			}
#endif	// VALIDATE

			Node1.m_MSA.Clear();
			Node2.m_MSA.Clear();

			delete[] Node1.m_Prof;
			delete[] Node2.m_Prof;
			Node1.m_Prof = 0;
			Node2.m_Prof = 0;
			}
		uTreeNodeIndex = GuideTree.NextDepthFirstNode(uTreeNodeIndex);
		}
	while (NULL_NEIGHBOR != uTreeNodeIndex);

	MakeRootMSA(Seqs, GuideTree, ProgNodes, Aln);
	delete[] ProgNodes;

#if	VALIDATE
	{
	unsigned uRootNodeIndex = GuideTree.GetRootNodeIndex();
	const ProgNode &RootProgNode = ProgNodes[uRootNodeIndex];
	AssertMSAEq(Aln, RootProgNode.m_MSA);
	}
#endif
	}
