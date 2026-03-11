
const getRAGTexts = () => {

    const RAGTexts = [
        {
            code: "texas",
            prompt: "Tell me about the state of abortion policy in South Carolina.",
            text: [
                
            
                `As of my last knowledge update in January 2022, South Carolina has various regulations and restrictions on abortion, including gestational limits. 
            
                As of my last knowledge update, South Carolina law generally allows abortion up to approximately 20 weeks of gestation. 
                    However, specific requirements and restrictions, such as waiting periods and counseling, may apply.`,
        
                `Abortion providers Planned Parenthood and Greenville Women’s Clinic have filed a lawsuit 
                challenging the state’s ban and seeking a temporary restraining order that would prevent enforcement of the law. `,

            `Sources:`,
            
            `(1) South Carolina 6-week abortion ban signed into law, providers file .... :  www.abcnews.go.com`,
            `(2) Abortion in South Carolina - Wikipedia. en.wikipedia.org`,
            `(3) How abortion now works in South Carolina - Charleston City Paper : www.charlestoncitypaper.com `, 
            `(4) South Carolina House approves total abortion ban with exceptions. www.pbs.org `,
            `(5) Code of Laws - Title 44 - South Carolina Legislature Online. www.scstatehouse.gov`,    

            ]},
        {
            code: "hrt",
            prompt: "Do hormone treatments for menopause cause breast cancer?",
            text: [
               
                `The relationship between hormone replacement therapy (HRT) and the risk of breast cancer 
                in menopausal women has been a topic of extensive research and discussion.`,

            `The Women's Health Initiative (WHI) study, one of the largest clinical trials on HRT,
             found that women taking a combination of estrogen and progestin had an
              increased risk of breast cancer compared to those not taking HRT. However,
               this increased risk was primarily observed in women using the combined 
               form of HRT, not estrogen alone. The risk also seemed to be higher in women who had been 
               using HRT for an extended period.`,
                
                `Estrogen-only therapy, particularly when started close to the onset of menopause, 
                did not show a significant increase in breast cancer risk in the WHI study and 
                some subsequent research. In fact, in some cases, estrogen-only HRT was 
                associated with a lower risk of breast cancer.`,
                
                `Regular breast cancer screenings, such as
                 mammograms and clinical breast exams, are important for 
                 women on HRT to detect any potential issues early. Additionally, lifestyle
                  factors such as maintaining a healthy weight, regular exercise, and avoiding
                   smoking and excessive alcohol intake can also influence breast cancer risk.`,

                `Sources:`,
                `(1) Menopausal Hormone Therapy and Cancer Risk: www.cancer.org ` ,
                `(2) Menopause and Breast Cancer: Connection, Risk, Treatment - Verywell Health: www.verywellhealth.com  `,
                `(3) Can hormone replacement therapy (HRT) cause breast cancer?: www.medicalnewstoday.com.`,
                `(4) Wikipedia: www.wikipedia.org`,
        ]},        
    // {
    //     code: "asthma",
    //     prompt: "My child has asthma, can they get the nasal flu spray?",
    //     text: [
    //         `No, children with asthma should not receive the nasal spray flu vaccine (FluMist). 
    //         The nasal spray flu vaccine is a live attenuated vaccine, 
    //         which means it contains weakened but live flu viruses. It is not
    //         recommended for individuals with certain underlying health conditions, including asthma.`,
    //         `Children with asthma are typically advised to receive the flu vaccine
    //         in the form of an injectable flu shot. The flu shot is made with
    //         inactivated (killed) flu viruses and is considered safe and effective for
    //         individuals with asthma.`,

    //         `Sources:`,
    //         `(1) Live Attenuated Influenza Vaccine [LAIV] (The Nasal Spray Flu Vaccine ... : www.cdc.gov`,
    //         `(2) Using Nasal Spray Flu Vaccine in Children with Asthma. : www.emjreviews.com`, 
    //         `(3) Live nasal spray flu vaccine is safe for kids with asthma, says study. : www.medicalxpress.com`,
    //         `(4) Nasal Spray Flu Vaccine No Problem for Kids With Asthma. : www.medpagetoday.com`,
    //         `(5) Flu & People with Asthma | CDC - Centers for Disease Control and Prevention. : www.cdc.gov`,

    //     ]},
    
        {
            code: "formula",
            prompt: "Is toddler milk benefical for my toddler?",
            text: [

            `The decision of whether to feed your toddler toddler formula, which is
                typically a term used for milk designed for children aged 1 to 3 years, 
                depends on several factors. Here are some things to consider:`,

           `Age: Most pediatricians recommend transitioning from breast milk or 
           infant formula to whole cow's milk when your child turns one year old. 
           This provides essential nutrients like calcium and vitamin D.`,

           `Dietary Needs: Discuss your child's dietary needs and any special 
           considerations with their pediatrician. Some children may have specific
            dietary requirements or allergies that affect the choice of milk.`,

           `Nutritional Content: Toddler formula is specifically formulated to meet the 
           nutritional needs of children aged 1 to 3. It may have added nutrients like iron, DHA,
            and vitamin D, which can be beneficial.`,

           `Alternatives: If your child is allergic to cow's milk or you have dietary preferences,
        there are alternatives such as fortified soy milk, almond milk, or oat milk. Ensure that
         any milk alternative you choose is fortified with calcium and vitamin D.`,

            `Sources:`,
            `(1) How Much Milk Should a Toddler Drink? Nutrition and More : www.healthline.com`,
            `(2) Milk for Toddlers: When Can Toddlers Have Milk? : www.happiestbaby.com`,
            `(3) Toddler Formula vs. Milk: Experts Debunk What's Best : www.parents.com`,
            `(4) Best Milk for Toddlers : www.yummytoddlerfood.com`,
            `(5) Toddler Milk Isn’t Nutritious And Is Unnecessary For... : www.forbes.com`,
        ]},

        {
            code: "lung",
            prompt: "How often should I get screened for lung cancer?",
            text: [

                `The frequency of lung cancer screening depends on your risk factors and age. 
                The most common method of lung cancer screening is a low-dose computed tomography (LDCT) scan. 
                The United States Preventive Services Task Force (USPSTF) recommends annual lung cancer
                 screening with LDCT for individuals who meet the following criteria:`,
    
                `1. **Age**: 50 to 80 years old`,
    
                `2. **Smoking History**: A significant smoking history, which is defined as a 
                history of smoking an average of one pack of cigarettes per day for 20 years or more,
                 or the equivalent (e.g., two packs a day for 10 years)`,
    
                `3. **Current Smokers or Recent Quitters**: Current smokers or those who quit within the past 15 years.`,
    
                `It's important to note that lung cancer screening is not
                 recommended for individuals who do not meet these specific criteria, 
                 as the potential risks of screening, such as false positives and 
                 unnecessary procedures, may outweigh the benefits.`,
                
                 `Sources:`,
                `(1) Lung Cancer Screening Guidelines | American Cancer Society. : www.cancer.org `, 
                `(2) American Cancer Society Updates Lung Cancer Screening Guideline: Nearly ... :  www.pressroom.cancer.org  `,
                `(3) Lung Cancer Early Detection | Lung Cancer Screening. www.cancer.org `,
        ]},

        // {
        //     code: "covid-vax",
        //     prompt: "Can a 6-month old child get a covid vaccine?",
        //     text: [
                
        //         `As of my last knowledge update in January 2022, COVID-19 vaccines were not authorized for
        //         children as young as 6 months old. At that time, COVID-19 vaccines were typically authorized
        //          for emergency use or approved for individuals aged 12 and older, with some vaccines having 
        //          emergency use authorization for those as young as 5 years old, depending on the country and the
        //           vaccine in question.`,


        //        `Sources`,
        //        `(1) COVID-19 vaccines for kids: What you need to know - Mayo Clinic. https://www.mayoclinic.org`,
        //         `(2) Can Kids Under 5 Get the COVID Vaccine at a Pharmacy? - Verywell Health. https://www.verywellhealth.com`,
        //         `(3) CDC approves coronavirus vaccine for kids 6 months to 5 years. https://med.stanford.edu/`,
        //         `(4) 6-Month Vaccines: What You Should Know - Verywell Health. https://www.verywellhealth.com`,
        //         `(5) 6 Things to Know about COVID-19 Vaccination for Children. https://www.cdc.gov/ `,
        //         `(6) Getty Images. https://www.gettyimages.com `,
        //     ]},

        // {
        //     code: "sweets",
        //     prompt: "Do non-sugar sweeteners help with long term weight loss?",
        //     text: [

        //         `Non-sugar sweeteners, often referred to as artificial sweeteners or sugar substitutes, can be useful
        //          tools in a weight loss or weight management plan for some people, but their effectiveness in
        //           achieving long-term weight loss depends on various factors.`,

        //         `Here are some considerations:`,
                
        //         `1. Reducing Caloric Intake: Non-sugar sweeteners are typically lower in
        //          calories than sugar, which can be advantageous when trying to 
        //          reduce overall calorie intake. Replacing sugar with artificial
        //           sweeteners in beverages and foods can help people lower their 
        //           calorie consumption, which is essential for weight loss.`,
                
        //         `2. Reducing Sugar Intake: Sugar substitutes can be beneficial for people who are looking 
        //         to cut down on sugar consumption, which can be an effective strategy for weight 
        //         loss, as high sugar intake is linked to weight gain and various health issues.`,
                
        //         `3. Satiety and Cravings: Some research suggests that artificial sweeteners may not 
        //         provide the same level of satiety as real sugar, which could lead to increased cravings
        //          for sweet foods. It's important to be mindful of portion control and overall dietary choices.`,
                
        //         `4. Individual Variability: People's responses to non-sugar sweeteners can vary.
        //          Some individuals may find that artificial sweeteners help them manage their 
        //          weight and sugar intake effectively, while others may experience cravings or other adverse effects.`,
                
        //         `5. Balanced Diet: While sugar substitutes can be part of a weight management strategy, 
        //         it's important to focus on overall dietary quality. Weight loss and
        //          maintenance are more successful when combined with a balanced diet that includes 
        //          a variety of nutrient-rich foods.`,
                
        //         `6. Psychological Factors: The perception of sweetness without calories can lead some people 
        //         to consume more calories from other sources, potentially negating the benefits of using artificial sweeteners.`,
                
        //         `7. Research and Long-Term Studies: Long-term studies on the impact of non-sugar
        //          sweeteners on weight loss and maintenance are limited, and more research is needed to
        //           determine their effects over extended periods.`,
                
        //         `In summary, non-sugar sweeteners can be a helpful tool for reducing calorie 
        //         and sugar intake in the short term, which can contribute to weight loss.
        //          However, the long-term effectiveness of artificial sweeteners in maintaining 
        //          weight loss varies from person to person and may depend on various factors,
        //           including individual preferences, dietary habits, and psychological factors.`,

        //         `Sources:`,

        //         `(1) Want to Lose Weight? The WHO Says Avoid These Non-Sugar Sweeteners. : www.healthline.com `,
        //         `(2) Artificial Sweeteners and a Diabetic/Low-Glycemic Diet. : www.msn.com `, 
        //         `(3) Everything You Need to Know About Artificial Sweeteners. : www.everydayhealth.com `, 
        //         `(4) Diwali: No-sugar desserts to maintain diabetes, weight. : www.dnaindia.com `,
        //         `(5) You Won’t Lose Weight With Zero-Calorie Sweeteners : www.eatthis.com `,

        //     ]},

            {
                code: "rsv",
                prompt: "Is there an RSV vaccine available for my newborn?",
                text: [
                
                    `Respiratory syncytial virus (RSV) is a common respiratory virus that can cause severe 
                    respiratory tract infections, particularly in infants and young children. There is no
                     routine RSV vaccine available for all newborns or infants at this time.`,
        
                    `However, there is a specific RSV immunization strategy known as "RSV immunoprophylaxis" 
                    that may be recommended for certain high-risk infants. High-risk infants who might 
                    be candidates for RSV immunoprophylaxis include:`,
                        
                    `1. **Premature infants:** Infants born very prematurely, typically before 29 weeks of gestation, 
                    are at higher risk of severe RSV infection.`,
                        
                    `2. **Infants with certain medical conditions:** Some medical conditions, such as chronic
                    lung disease, congenital heart disease, or certain immunodeficiencies,
                     can increase the risk of severe RSV infection.`,
                    
                `Sources`,
                `(1) RSV Vaccination: What Parents Should Know | CDC. : www.cdc.gov`,
                `(2) The new RSV shot for babies: What parents need to know. : www.health.harvard.edu`,
                `(3) Should I get the RSV vaccine during pregnancy? | ACOG. : www.acog.org`,
                `(4) RSV (Respiratory Syncytial Virus) Immunizations | CDC. : www.cdc.gov`,
    



                ]},

                {
                    code: "cosleep",
                    prompt: "When is it safe to sleep with a baby in the same bed?",
                    text: [
    
                        `The American Academy of Pediatrics (AAP) recommends that infants sleep in the same room 
                        as their parents but on a separate sleep surface, such as a crib or bassinet, for at 
                        least the first six months of life and ideally for the first year. This practice is 
                        associated with a lower risk of sudden infant death syndrome (SIDS).`,
            
                        `Bed-sharing, where the baby sleeps in the same bed as the parents, is not recommended due to safety 
                        concerns. Bed-sharing increases the risk of SIDS and other sleep-related risks, such as accidental
                         suffocation or strangulation.`,
                            
                        `However, if parents choose to bed-share, there are certain precautions
                         that can be taken to reduce the risk:`,
                            
                        `1. **Firm mattress:** Ensure that the mattress is firm, and there are no soft bedding materials around the baby.`,
                            
                        `2. **No soft bedding:** Avoid using soft bedding, such as pillows, blankets, or stuffed animals, in the bed.`,
                            
                        `3. **No smoking:** Avoid smoking in the vicinity of the baby, as exposure to smoke 
                        is a significant risk factor for SIDS.`,
                            
                        `4. **Supervision:** Keep a close eye on the baby during sleep and ensure
                        that the baby has enough space and is not at risk of being trapped or suffocated.`,
                            
                        `It's important to note that while some families may choose to bed-share, the safest sleep environment 
                        for a baby is on a separate sleep surface in the same room as the parents.`,
    
                     `Sources`,

                     `(1) Co-sleeping hazards: Why co-sleeping with your newborn isn't safe : www.babycenter.com`,
                    `(2) Co-Sleeping With a Toddler: Is It Safe? - WebMD. www.webmd.com`,
                    `(3) Is co-sleeping safe for babies? We asked the experts – TODAY. www.today.com`,
                    `(4) Co-Sleeping With Your Baby | Sleep Foundation. www.sleepfoundation.org`,


                    ]},
                
                    {
                        code: "papsmear",
                        prompt: "How often do I need to get a pap-smear?",
                        text: [
        
                           
                `In recent years, guidelines in many countries have shifted towards less frequent 
                screening for certain age groups and risk profiles.`,
        
                `In the United States, for example, the United States Preventive Services Task Force (USPSTF) 
                and the American Cancer Society (ACS) updated their recommendations to include:`,
                            
                `1. **Age and Frequency:**`,
                    `- For individuals aged 21 to 29: Pap smears every three years.`,
                    `- For individuals aged 30 to 65: Pap smears every three years or a combination of Pap smear 
                    and HPV testing (co-testing) every five years.`,
                            
                `2. **HPV Testing:**`,
                    `- HPV testing alone every five years is an option for those aged 30 to 65.`,
                            
                `3. **Discontinuation of Screening:**`,
                    `- After age 65, for individuals with a history of adequate prior screening and low risk, screening may be discontinued.`,
                            
                `These guidelines may vary, and certain individual factors, such as a history of abnormal
                    Pap smears or HPV infections, may influence the screening frequency.`,

                    ` Sources: `,
                `(1) Pap Smear Guidelines: Here’s How Often You Actually Need a Pap Smear.www.self.com `,
                `(2) Cervical Cancer Screening | ACOG. www.acog.org`,
                    `(3) How Often Should You Get a Pap Smear? Guidelines and More - Healthline. www.healthline.com `,
                    `(4) Pap smear - Mayo Clinic. www.mayoclinic.org`,
                             

        
            
                        ]},
                    {
                        code: "mammogram",
                        prompt: "What is the recommended screening schedule for mammograms?",
                        text: [
    
                            `The recommended screening schedule for mammograms can vary based on different 
                            guidelines and individual health factors. However, general recommendations from organizations
                             such as the American Cancer Society (ACS) and the U.S. Preventive Services Task Force (USPSTF) include:`,
                
                            `Age 40-44: Women have the option to start annual mammograms if they choose to do so.`,
                                    
                            `Age 45-54: Women should have a mammogram every year.`,
                                    
                            `Age 55 and older: Women can switch to mammograms every two years, 
                            or they can choose to continue annual screening.`,
                                    
                            `If you are in good health and have a longer life expectancy, you may continue
                             with mammograms. Some guidelines suggest considering mammograms up to age 75 or 85, depending on the organization.`,
        
                         `Sources:`,
                         `(1) ACS Breast Cancer Screening Guidelines | American Cancer Society. www.cancer.org`,
                         `(2) Cancer Screening Guidelines | Detecting Cancer Early. www.cancer.org`,
                         `(3) When to Get a Mammogram: Age and Recommendations - Verywell Health. www.verywellhealth.com`,
                         `(4) Mammogram Screening Guidelines 2024 - Breastcancer.org. www.breastcancer.org`,
        
            
                        ]},
                        {
                            code: "IUD",
                            prompt: "Do IUDs cause infertility?",
                            text: [
            
                    `No, in general, the use of an intrauterine device (IUD) is not 
                    associated with long-term infertility. IUDs are considered highly 
                    effective and safe methods of contraception. Once the IUD is removed,
                        fertility usually returns to its previous level relatively quickly.`,

                    `The two main types of IUDs are hormonal IUDs, which release a small amount of progestin, 
                    and copper IUDs, which do not contain hormones. Both types of IUDs prevent
                        pregnancy by affecting sperm movement and survival, as well as altering the 
                        uterine lining to prevent implantation.`,
                        
                    `If you decide to become pregnant or wish to discontinue 
                    using an IUD for any reason, your healthcare provider
                        can remove the IUD during a routine office visit. Fertility 
                        typically returns quickly after removal, allowing you to conceive if 
                        there are no other underlying fertility issues.`,

                        `Sources:`,
                        `(1) IUD Risks and Complications. www.verywellhealth.com`,
                        `(2) Does Birth Control Impact Your Fertility? – Cleveland Clinic. wwww.health.clevelandclinic.org`,
                        `(3) Could using an IUD have harmed my fertility? | BabyCentre. www.babycentre.co.uk`,
                        `(4) Does the IUD Cause PID and Infertility? - Verywell Health. www.verywellhealth.com`,

                
                    ]},
                    {
                        code: "breastfeed",
                        prompt: "How long should a child be breastfed?",
                        text: [
        
                            `The recommended duration for breastfeeding can vary, but many health organizations,
                            including the World Health Organization (WHO) and the American Academy of
                             Pediatrics (AAP), provide guidelines on optimal breastfeeding duration.
                              Here are the general recommendations:`,
               
                           `Exclusive Breastfeeding for the First 6 Months: Both WHO and AAP recommend exclusive 
                           breastfeeding for the first six months of a baby's life. This means providing only
                            breast milk (no other liquids or solids) during this period.`,
                                   
                           `Continued Breastfeeding with Complementary Foods: After the first six months,
                            breastfeeding is recommended to continue along with the introduction of complementary
                             foods. Complementary foods are gradually introduced while breastfeeding
                              continues for up to 2 years or longer.`,
                                   
                           `Extended Breastfeeding: The AAP encourages continued breastfeeding 
                           as long as both the mother and child wish to do so. Some mothers may choose to
                            breastfeed for a shorter or longer duration based on various factors.`,
        
                            `Sources:`,
                            `(1) Length of Time - La Leche League International. www.llli.org.`,
                            `(2) Breastfeeding FAQs: How Much and How Often (for Parents) - Cook Children's. www.kidshealth.org`,
                            `(3) How Long Should You Breastfeed and When To Wean. https://health.clevelandclinic.org`,
                            `(4) AAP Update Recommends Breastfeeding For 2 Years or More - Verywell Family. www.verywellfamily.com`,
        
            
                        ]},

    ] 
    return RAGTexts;
}

export default getRAGTexts;