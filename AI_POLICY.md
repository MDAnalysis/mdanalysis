# 🤖 MDAnalysis AI-generated contributions policy 🤖

* Version: **2.0 (2026-11-30)**
* License: **[CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/)**

This document defines the [MDAnalysis organization's](https://github.com/MDAnalysis/) policy regarding AI-generated content. This policy applies to all aspects of the MDAnalysis project, including all the GitHub repositories under the organization unless a project defines its own explicit AI Policy.

For the purposes of this document, _content_ is defined as any contribution to the MDAnalysis organization, both code and non-code, including but not limited to:
- Code submitted for inclusion in MDAnalysis, or other organization-owned repositories
- Documentation of any kind
- Communication content (e.g. issues, discussions, pull request descriptions, and reviews)

We recognize that the use of AI in open source software raises many difficult questions and acknowledge that views, approaches, and technical capabilities are rapidly changing. 
We invite comments and further discussion (via our [discussion channels](https://www.mdanalysis.org/pages/community/)) and expect to periodically review this policy.


## Guiding principles

* Humans are always responsible for AI output.
* Contributors must be transparent and honest; this includes that AI-use must always be declared.
* In forums designed for humans (e.g., issue tracker, discussions) only humans are allowed to speak. We want to spend our time talking with humans, not with LLMs (or human proxies for LLMs).
* Human time is precious. Core Developers can always decide that a contribution does not justify the review effort and close it.
* The project benefits most from expertise, not just code/documentation.
* Scientific code has to be correct. If no human can assess or validate the code then we cannot include it.
* Best practice norms, e.g. software contribution standards set by the project, must always be followed.


## Policy foundations

MDAnalysis allows AI generated content in code and documentation but not in forums intended for humans (with some limited exceptions). 

**Any form of AI used must always be declared**.

### Expertise and trust matter

The project benefits from human expertise but not from low-effort contributions and hence CoreDevelopers will exercise broad discretion in rejecting PRs outright that appear to be AI-generated slop or situations where the AI generated code is not fully understood by the contributing author. This may mean that PRs are treated differently, depending on who is submitting them. Someone who has contributed to the project over time may be trusted to use generative AI responsibly and have their PR reviewed while a new contributor’s PR is just rejected.

We want to make explicit that trust factors into our decisions. Trust is something that we are building together through repeated, positive interactions in a range of modalities and we do not consider it as something one can shortcut or outsource to an algorithm.


### How contributions will be assessed

Please understand that reviewing capacity is very limited for the project, so large PRs which appear to not have the requisite understanding might not get reviewed, and eventually closed or redirected.
Understanding from authors about what code they are contributing and how it fits into the larger picture is crucial; for new features in particular, we are looking for domain expertise.
If reviewers feel there is a lack of understanding from looking at the code, PR descriptions, or just general interactions with a contributor, we reserve the right to close the PR to save time and resources.

Given the limited reviewing capacity, time and effort of reviewing and educating contributors will be more targeted to PRs of initial high quality and/or contributors where interactions are constructive and of good quality.
This means that for new or infrequent contributors, we may refuse reviews of PRs without an associated initial discussion on a corresponding issue. We look for a high-level description and motivation of what a PR will include. We do not recommend opening PRs without having had a discussion with a member of the [Core Developer team](https://www.mdanalysis.org/pages/governance/#mdanalysis-core-developers) on an associated issue.

## Policy details

### 1. AI usage must always be declared

All contributions must be fully transparent about any AI usage. Failure to disclose AI use (e.g., in pull request description details)
will be considered non-compliance with this policy.

If we suspect that AI tools were used as part of a contribution, MDAnalysis developers may ask you for clarification.
Please do not take offense if this happens to you, distinguishing between AI-generated and human-authored content is increasingly difficult!

Each PR that used generative AI tool(s) should contain an "assisted-by" tag in the commit message that indicates the model and version, 
```
Assisted-by: AGENT_NAME:MODEL_VERSION [tool1] [tool2]
```
for instance:
```
Assisted-by: ChatGPT:GPT-6.1-Sol-Light

Assisted-by: Claude:claude-opus-5-5 Cursor
```

(Ask the model to generate the appropriate tag for you or provide your best guess.)


### 2. Restricted use of AI in forums

A "forum" is a place built for human communication; this includes (but is not limited to) issue trackers, PR comments, Discord, GitHub discussions. We require humans to write in these forums and do not want to engage with generated content.

You may use AI help to, for instance, generate a reproducer case, or format tables to communicate results. But a human has to write the overall post and be in charge of the content.



### 3. Language translation

MDAnalysis aims to foster an inclusive environment for all contributors and users. We recognise that language barrier can limit participation
in open-source projects, particularly for non-native English speakers.

Whilst we prefer you communicate with us in your own words (and languages if needed), we acknowledge that AI-based translation tools may sometimes
be necessary for communicating in issues / discussions. If such tools are used, we ask you to disclose their use and purpose.

### 4. Contributors are responsible authors and must be aware of, and declare, any derivative work

Contributors are responsible for the content and provenance of every aspect of their contribution, including any design decisions.
Contributors must be aware of whether their work derives from external sources / codebases and must declare any derivative work.

Code reviewers may ask you to explicitly confirm that your code does not infringe on any copyright and licenses and in this case the burden is on you to address this question.

### 5. Human reviewers are required

All code merged into MDAnalysis repositories must be reviewed by a human reviewer. Instructions / suggestions from human reviewers always take precedence over those of non-human reviewers.

### 6. Consequences for failure to comply to this policy

Any failures to comply with this policy should be reported to the MDAnalysis core developers through
our [standard communications channels](https://www.mdanalysis.org/pages/community/#ask-questions--get-help) or email to mdanalysis@numfocus.org .
The Core Developers may enforce this policy through a range of actions, which may include restricting or banning individuals from future contributions.

There is a zero tolerance policy for intentionally failing to disclose or misrepresent the extent of AI usage.

## Acknowledgements

We acknowledge the [Kornia AI policy](https://github.com/kornia/kornia/blob/main/AI_POLICY.md) and discussions
within the [Scientific Python Community](https://discuss.scientific-python.org/t/a-policy-on-generative-ai-assisted-contributions/1702/18)
for providing guidance and insight in shaping this policy.

