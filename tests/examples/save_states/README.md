The origin of the files in this folder:

- `amplify_4_all_example_primers.tsv` is from
  `Amplify4/Amplify4.app/Contents/Resources/Amplify Sample Primers.primers`
- `amplify_4_example_template.txt` is from
  `Amplify4/Amplify4.app/Contents/Resources/Amplify Sample Target.rtf`
- `amplify_4_full_example.yaml` was created by loading in
  `amplify_4_all_example_primers.tsv` and `amplify_4_example_template.txt` and
  saving the state.
- `simple.yaml` was created from Amplify 4's samples.
- `long_amplicons.yaml` was created from Amplify 4's samples, it contains a
  template and primer set that produce long amplicons with high raw $Q$ scores.
- `lots_of_amplicon.yaml` was created from Amplify 4's samples contains a
  template and primer set that produce many amplicons across a wide range of $Q$
  scores.
- `overlapping_primers.yaml` was supplied by Petr Vasilchikov, and it contains
  primers which are overlapping.
