# Generated by Django 3.0.5 on 2020-04-03 14:00

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('drugstone', '0022_auto_20200403_1319'),
    ]

    operations = [
        migrations.AddField(
            model_name='protein',
            name='proteins',
            field=models.ManyToManyField(related_name='_protein_proteins_+', through='drugstone.PPI', to='drugstone.Protein'),
        ),
        migrations.AlterUniqueTogether(
            name='ppi',
            unique_together={('from_protein', 'to_protein')},
        ),
    ]
